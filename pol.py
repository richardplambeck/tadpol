# ---------------------------------------------------------------------------------------------------------- #
     # note that phiDeg=0 corresponds to no perp to plane of incidence, apparently
# pol.py

#from Numeric import *
import numpy 
import math
import cmath
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

# ---------------------------------------------------------------------------------------------------------- #
# define basis vectors in reference coordinate system (aligned with OMT axes)
# each vector consists of 2 complex numbers [ (Re_x,Im_x), (Re_y,Im_y) ]
# R and L propagate in +Z direction according to right hand rule
# note: x.conjugate() gives complex conjugate of a number
# ---------------------------------------------------------------------------------------------------------- #
X = numpy.array( [1,0], dtype=complex )
Y = numpy.array( [0,1], dtype=complex )
R = numpy.array( [1./math.sqrt(2.), -1j*(1./math.sqrt(2.))] ) 
L = numpy.array( [1./math.sqrt(2.), +1j*(1./math.sqrt(2.))] ) 
clight = 29.9792458	# speed of light, cm/nanosec

# ---------------------------------------------------------------------------------------------------------- #
# returns new basis vector in a coordinate system rotated by thetaDegrees
# ---------------------------------------------------------------------------------------------------------- #
def Jrot ( vec, thetaDegrees ) :
  rad = math.pi * thetaDegrees/180.
  rotmat = numpy.array( [[math.cos(rad),math.sin(rad)],[-math.sin(rad),math.cos(rad)]] )
  return numpy.dot(rotmat,vec)

# ---------------------------------------------------------------------------------------------------------- #
# returns basis vector after passing through polarizer section
# component 2 (Y-axis) is advanced by delayDegrees relative to component 1 (X-axis)
# ---------------------------------------------------------------------------------------------------------- #
def Jdelay ( vec, delayDegrees ) :
  rad = math.pi * delayDegrees/180.
  rotmat = numpy.array( [ [1, 0], [0, cmath.exp(-1.j * rad)] ] )
  return numpy.dot( rotmat,vec )

# --- delays both components --- #
def J2delay ( vec, delayDegrees ) :
  rad = math.pi * delayDegrees/180.
  rotmat = numpy.array( [ [cmath.exp(-1.j * rad), 0], [0, cmath.exp(-1.j * rad)] ] )
  return numpy.dot( rotmat,vec )

# --- delays by different amounts --- #
def J3delay (vec, delay1Deg, delay2Deg ) :
  rad1 = math.pi * delay1Deg/180.
  rad2 = math.pi * delay2Deg/180.
  rotmat = numpy.array( [ [cmath.exp(-1.j * rad1), 0], [0, cmath.exp(-1.j * rad2)] ] )
  return numpy.dot( rotmat,vec )

# ---------------------------------------------------------------------------------------------------------- #
# returns basis vector after transmission through a beamsplitter
# tpel is thickness in inches, angI is angle of incidence in degrees
# X axis is assumed parallel to the plane of incidence; light polarized parallel to the plane of incidence
#   is less strongly reflected; at Brewster's angle none will be reflected
# Y axis is perp to plane of incidence; it is more strongly reflected
# ---------------------------------------------------------------------------------------------------------- #
def Jbsplit ( vec, tpel=.001, angI=45., fGHz=230. ) :
   [tpar,tperp,Rpar,Rperp] = pellicle( tpel=tpel, angI=angI, fGHz=fGHz)
   rotmat = numpy.array( [ [tpar, 0.], [0., tperp] ] )
   return numpy.dot( rotmat,vec ) 

# ---------------------------------------------------------------------------------------------------------- #
# returns cutoff freqs [fcX, fcY], in GHz, for faceted circular waveguide that is squeezed along Y-direction
# r is waveguide radius in inches, x is the normalized facet depth f/r
# based on 5th order polynomial fit to analytic results from Appendix A, or HFSS simulations
# use source='HFSS' for HFSS (default), source='Wang' for analytic fits
# ---------------------------------------------------------------------------------------------------------- #
def cutoff( r, x, source='HFSS' ) :
  if ( source == 'Wang') :
    kX = 1.841184 + 0.309847*x + 8.60483*pow(x,2.) - 29.3916*pow(x,3.) + 75.4104*pow(x,4.) - 67.4144*pow(x,5.)
    kY = 1.841184 - 0.0900839*x - 3.30519*pow(x,2.) + 13.3106*pow(x,3.) - 26.8144*pow(x,4.) + 22.970*pow(x,5.)
  else :
    kX = 1.841184 + 0.301574*x + 8.9118*pow(x,2.) - 33.253*pow(x,3.) + 93.2359*pow(x,4.) - 94.615*pow(x,5.)
    kY = 1.841184 - 0.0862305*x - 3.41638*pow(x,2.) + 14.65*pow(x,3.) - 32.7615*pow(x,4.) + 31.7498*pow(x,5.)
  fcX = clight * kX / (2. * math.pi * r * 2.54)
  fcY = clight * kY / (2. * math.pi * r * 2.54)
  return [fcX, fcY]
     
# ---------------------------------------------------------------------------------------------------------- #
# compute length in inches of faceted guide to achieve differential phase shift dphi0, in degrees 
# r is waveguide radius in inches, f is facet depth in inches, fGHz is freq in GHz
# fcX and fcY are override cutoff freqs in GHz; if either is zero, cutoff freqs will be taken from 
#   subroutine cutoff, using 5th order polynomial fit to the HFSS-derived cutoffs
# ---------------------------------------------------------------------------------------------------------- #
def length( r=0.0235, f=0.006, dphi0=90., fGHz=230., fcX=0., fcY=0. ) :
  if (f == 0.) : return 0.
  if ( (fcX == 0.) or (fcY == 0.) ) :
    [fcX,fcY] = cutoff( r, f/r, source='HFSS' )
  length = clight * dphi0 / (2.54 * 360.* (math.sqrt(fGHz*fGHz-fcY*fcY) - math.sqrt(fGHz*fGHz-fcX*fcX)))
  return length

# --- L90 retarder lengths in table 1 are from the actual HFSS freqs, not the polynomial fit --- #
def table1 () :
  print "0.001 %8.4f" % length( r=0.0235, fcX=149.307, fcY=146.471, fGHz=230.)
  print "0.002 %8.4f" % length( r=0.0235, fcX=153.119, fcY=145.203, fGHz=230.)
  print "0.003 %8.4f" % length( r=0.0235, fcX=158.048, fcY=143.674, fGHz=230.)
  print "0.004 %8.4f" % length( r=0.0235, fcX=163.985, fcY=142.030, fGHz=230.)
  print "0.005 %8.4f" % length( r=0.0235, fcX=170.938, fcY=140.363, fGHz=230.)
  print "0.006 %8.4f" % length( r=0.0235, fcX=178.985, fcY=138.732, fGHz=230.)
  print "0.007 %8.4f" % length( r=0.0235, fcX=188.256, fcY=137.176, fGHz=230.)

# --- generate smooth curves of cutoff freq for FIG 4 --- #
def fig4() :
  r = .0235
  for f in arange( 0.000, 0.0072, 0.0002) :
    [fcX,fcY] = cutoff( r, f/r )
    print "%6.4f  %8.4f  %8.4f" % (f, fcX, fcY)
   
# ---------------------------------------------------------------------------------------------------------- #
# compute phase delay in degrees of Y-pol vs X-pol signals traveling through length L; dimensions in inches
# fcX and fcY are override cutoff freqs in GHz; if either is zero, use values from subroutine cutoff
# ---------------------------------------------------------------------------------------------------------- #
def dphi( r=0.0235, f=0.006, L=0.001, fGHz=230., fcX=0., fcY=0. ) :
  if ( (fcX == 0.) or (fcY == 0.) ) :
    [fcX,fcY] = cutoff( r, f/r, source='HFSS' )
  dphi = 360. * L * 2.54 * (math.sqrt(fGHz*fGHz-fcY*fcY) - math.sqrt(fGHz*fGHz-fcX*fcX)) / clight
  return dphi

# ---------------------------------------------------------------------------------------------------------- #
# compute [length, differential phase shift] through a radiused adapter from faceted circ to circ guide
# Rc is the cutter radius in inches; axis of cutter is perpendicular to axis of waveguide
# ---------------------------------------------------------------------------------------------------------- #
def dphitaper2 ( r=0.0235, f=0.006, Rc=0.125, fGHz=230., nsteps=1000 ) :
  df = f/nsteps                                                  # increment in facet depth
  f1 = f                                                         # facet depth at start of segment
  x1 = 0.                                                        # x-coordinate at start of segment
  totphase = 0.
  for n in range(nsteps) :
    f2 = f1 - df                                                 # facet depth at end of segment
    x2 = sqrt(Rc*Rc - pow((f - f2 - Rc), 2.))                    # x-coordinate at end of segment
    dphase = dphi(r=r, f=(f1+f2)/2., L=(x2-x1), fGHz=fGHz )      # using avg facet depth
    totphase = totphase + dphase
    x1 = x2
    f1 = f2
  return [x2, totphase]           # return X length of transition in inches, and total diff phase through it

# ---------------------------------------------------------------------------------------------------------- #
# predictions for Kband scale models; gap = 0.382" for straight case 1, 0.338" for straight case 2
# ---------------------------------------------------------------------------------------------------------- #
def fig8ab () :
  ofile = open("KbandStraight.dat","w")
  for f in arange(19.,30.05,.05) :
    ofile.write("%6.2f %8.4f %8.4f\n" % 
     (f, dphi( r=0.2275, f=.0365, L=1.822, fGHz=f ), dphi( r=0.2275, f=.0585, L=1.822, fGHz=f )) )
  ofile.close()

def fig8c () :
  ofile = open("KbandCurved.dat","w")
  for f in numpy.arange(19.,30.05,.05) :
    [x2, ph] = dphitaper2( r=0.2275, f=.05675, Rc=1.210, fGHz=f )    # min gap is 0.3415"
    ofile.write("%6.2f %8.4f \n" % (f, 2.*ph))
  ofile.close()
  
# ---------------------------------------------------------------------------------------------------------- #
# reflection coefficient from circular guide of radius r1 (inches) to guide of radius r2 (inches)
# ---------------------------------------------------------------------------------------------------------- #
def reflect( fGHz, r1, r2 ) :
  fc1 = clight/(3.4126 * r1 * 2.54)
  fc2 = clight/(3.4126 * r2 * 2.54)
  Z1 = fGHz * 377./math.sqrt(fGHz*fGHz - fc1*fc1)
  Z2 = fGHz * 377./math.sqrt(fGHz*fGHz - fc2*fc2)
  vR = (Z1 - Z2)/(Z1 + Z2)
  print "voltage, power reflection coefficients: %.3e %.3e" % (vR, vR*vR)

# ---------------------------------------------------------------------------------------------------------- #
# check algebra of example at end of section 5
# ---------------------------------------------------------------------------------------------------------- #
def example() :
  v1 = X
  v2 = Jrot(v1, 45)
  v3 = Jdelay(v2, 90)
  v4 = Jrot(v3, -45)
  print v4
  rad = math.pi/4
  print cmath.exp(+1.j * rad) * v4
  print "D_L ", numpy.dot(v4,L), abs(numpy.dot(v4,L))	# note: R* = L
  print "D_R ", numpy.dot(v4,R), abs(numpy.dot(v4,R))	# note: L* = R

# ---------------------------------------------------------------------------------------------------------- #
# writes out dimensions of a single polarizer to 1 line of output file; 
# convert inches to mils for more compact format 
# ---------------------------------------------------------------------------------------------------------- #
def dumpDimensions( dimfile, r, a1, f1, L1, a2, f2, L2, a3, f3, L3 ) :
  ofile = open(dimfile, "a")
  ofile.write(" %6.3f %7.2f %5.3f %6.2f" % (1000.*r, a1, 1000.*f1, 1000.*L1) )
  ofile.write(" %7.2f %5.3f %6.2f" % (a2, 1000.*f2, 1000.*L2) )
  ofile.write(" %7.2f %5.3f %6.2f\n" % (a3, 1000.*f3, 1000.*L3) )
  ofile.close()

# ---------------------------------------------------------------------------------------------------------- #
# writes out one polarizer dimension (e.g., facet depth of section 1) per line for multiple polarizers 
# ---------------------------------------------------------------------------------------------------------- #
def dumplist ( ofile, x, label, fmt, mult ):
  ofile.write( label ) 
  ntrials = len(x)
  for n in range(ntrials) :
    ofile.write( fmt % (mult * x[n] ) )
  ofile.write("\n")

# ---------------------------------------------------------------------------------------------------------- #
# compute ampX, phsX, ampY, phsY, phsdif given a pol vector [ (re X, im X), (re Y, im Y) ]
# ---------------------------------------------------------------------------------------------------------- #
def ampPhs( vec ) :
  phsX = 180. * math.atan2(vec[0].imag,vec[0].real) / math.pi
  phsY = 180. * math.atan2(vec[1].imag,vec[1].real) / math.pi
  phsdif = phsX - phsY
  if (phsdif > 180.) :
    phsdif = phsdif - 360.
  if (phsdif < -180.) :
    phsdif = phsdif + 360.
  return [ abs(vec[0]), phsX, abs(vec[1]), phsY, phsdif ]

# ---------------------------------------------------------------------------------------------------------- #
# compute leakages vs freq for 1-section, 2-section, or 3-section polarizers
# read polarizer dimensions from 'dimfile', write leakages to 'leakfile'
# note: one polarizer per ROW on input, one per COLUMN on output
# dimfile should have 10 columns (as in dumpDimensions); angles in degrees, lengths in mils
# optional: apfile lists amps and phase difference of X and Y components for comparison with HFSS
# optional: evaluate leakage after beamsplitter  apel,tpel,aIpel,aeval
#    apel = angle of plane of incidence, relative to X=0 axis 
#    tpel = beamsplitter ("pellicle") thickness in inches
#    aIpel = angle of incidence to the beamsplitter, degrees (0 = beamsplitter normal to axis)
# optional: evaluate leakage at angle aeval (shouldn't affect magnitude of the leakage) 
# optional: compute correlation efficiency < V1 V0* > for all polarizers vs the first one
# ---------------------------------------------------------------------------------------------------------- #
def computeLeakage( dimfile, leakfile, apfile=None, apel=0., tpel=0., aIpel=0., aeval=0., efficfile=None ) :
  r  = []   # empty list of circular waveguide radii
  a1 = []   # angle of section 1 relative to input Y-pol
  f1 = []   # facet depth of section 1
  L1 = []   # length of section 1
  a2 = []   # angle of section 2 relative to section 1
  f2 = []   # facet depth of section2
  L2 = []   # length of section 2
  a3 = []   # angle of section 3 relative to section 2
  f3 = []   # facet depth of section 3
  L3 = []   # length of section 3
  
  # --- read the dimensions into internal lists --- #
  infile = open(dimfile, "r")
  ofile = open(leakfile, "w")
  if (efficfile) : ofile3 = open( efficfile, "w" )
  if (apfile) :
    ofile2 = open( apfile, "w")
    ofile2.write("# fGHz, abs(Ex), phs(Ex), abs(Ey), phs(Ey), phsdif, leakage, phsleak\n" )
  ntrials = 0
  for line in infile:
    if line.startswith("#"):
      ofile.write(line)
    else:
      a = line.split()                      # split line into string tokens
      r.append( float(a[0])/1000. )         # convert diameter from mils back to inches
      a1.append( float(a[1]) )
      f1.append( float(a[2])/1000. )	    # convert facet depth from mils back to inches
      L1.append( float(a[3])/1000. )
      a2.append( float(a[4]) )
      f2.append( float(a[5])/1000. )	    # convert facet depth from mils back to inches
      L2.append( float(a[6])/1000. )
      a3.append( float(a[7]) )
      f3.append( float(a[8])/1000. )	    # convert facet depth from mils back to inches
      L3.append( float(a[9])/1000. )
      ntrials = ntrials + 1
  infile.close()
    
  # --- compute leakage array for each set of polarizer dimensions --- #
  freq = numpy.arange(200.,271.)
  leakage = numpy.zeros([len(freq),ntrials],Float)    # create array to hold leakages
  effic = numpy.zeros([len(freq),ntrials],Float)      # create array to hold correlation efficiencies
  for n in range(ntrials) :
    m = 0                                       # freq index
    for fGHz in freq :

      v1 = Y                                    # Y-pol incident on section 1
      v2 = Jrot( v1, a1[n] )
      afinal = -1.*a1[n]
      v3 = Jdelay( v2, dphi( r[n], f1[n], L1[n], fGHz ) )

      if (f2[n] > 0.) :                         # section 2, if present
        v2 = Jrot( v3, a2[n] )
        afinal = afinal - a2[n]
        v3 = Jdelay( v2, dphi( r[n], f2[n], L2[n], fGHz ) )

      if (f3[n] > 0.) :                         # section 3, if present
        v2 = Jrot( v3, a3[n] )
        afinal = afinal - a3[n]
        v3 = Jdelay( v2, dphi( r[n], f3[n], L3[n], fGHz ) )

      v1 = Jrot(v3, afinal)	             	# rotate back to original reference frame

      if (tpel > 0.) :		                # optional beamsplitter section
        v2 = Jrot( v1, apel )                   # rotate X axis parallel to plane of plane of incidence
        v3 = Jbsplit(v2, tpel, aIpel, fGHz )    # apply amplitude and phase shifts
        v1 = Jrot( v3, -1.*apel )		# rotate back to original reference frame

      vout = Jrot( v1, aeval )			# evaluate leakage at angle aeval

      if (n == 0) : vsave = numpy.array( [vout[0].conjugate(), vout[1].conjugate() ] )
        # save polarization vector of 1st trial for optional efficiency calculations

      # --- assume RCP out (true for positive angles a1,a2,a3) --- #
      cleak = numpy.dot(vout, R)                      # single complex number; note R = L*
      leakage[m,n] = abs(cleak)	                # magnitude of the leakage 

      if (efficfile) : effic[m,n] = abs(numpy.dot(vsave,vout))
      if (apfile) :
        [ampX, phsX, ampY, phsY, phdif ] = ampPhs( vout )
        ofile2.write("%6.1f  %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (fGHz, ampX, phsX, ampY, phsY, leakage[m,n]))
      m = m + 1

  if (apfile) : ofile2.close()

  # --- write dimensions in columns to output file --- #
  dumplist ( ofile, r,  "# r  :", " %6.3f", 1000.) 
  dumplist ( ofile, a1, "# a1 :", " %6.2f", 1.) 
  dumplist ( ofile, f1, "# f1 :", " %6.3f", 1000.) 
  dumplist ( ofile, L1, "# L1 :", " %6.2f", 1000.) 
  dumplist ( ofile, a2, "# a2 :", " %6.2f", 1.) 
  dumplist ( ofile, f2, "# f2 :", " %6.3f", 1000.) 
  dumplist ( ofile, L2, "# L2 :", " %6.2f", 1000.) 
  dumplist ( ofile, a3, "# a3 :", " %6.2f", 1.) 
  dumplist ( ofile, f3, "# f3 :", " %6.3f", 1000.) 
  dumplist ( ofile, L3, "# L3 :", " %6.2f", 1000.) 
  if (tpel == 0.) :
    ofile.write("# evaluated without any beamsplitter\n")
  else :
    ofile.write("# beamsplitter at angle apel = %.2f\n" % apel)
    ofile.write("# beamsplitter thickness tpel = %.4f\n" % tpel)
    ofile.write("# beamsplitter angIncidence aIpel = %.2f\n" % aIpel)
  ofile.write("# leakages evaluated after rotation by aeval= %.2f\n" % aeval)

  # --- dump leakages at each freq, for all polarizers; also compute avg and rms --- #
  m = 0
  for fGHz in freq :
    ofile.write("%6.1f" % fGHz)
    if (efficfile) : ofile3.write("%6.1f" % fGHz)
    sum = 0.
    for n in range(ntrials) :
      ofile.write(" %6.4f" % leakage[m,n])
      sum = sum + leakage[m,n]
      if (efficfile) : ofile3.write(" %6.4f" % effic[m,n] )
    avg = sum/ntrials
    var = 0.
    for n in range(ntrials) :
      dif = leakage[m,n] - avg
      var = var + dif*dif
    if (ntrials > 1) :
      var = var/(ntrials-1)
    ofile.write("   %6.4f %6.4f\n" % (avg,math.sqrt(var)))
    if (efficfile) : ofile3.write("\n")
    m = m + 1
  ofile.close()     
  if (efficfile) : ofile3.close()

# ---------------------------------------------------------------------------------------------------------- #
# leakage of idealized polarizers
# ---------------------------------------------------------------------------------------------------------- #
def fig9() :
  r = .0235
  f = 0.006
  L90 = length( r=r, f=f, dphi0=90., fGHz=230. ) 
  L180 = 2. * L90

  a1 = 45.
  dumpDimensions( 'dim9A.dat', r, a1, f, L90, 0., 0., 0., 0., 0., 0.)
  computeLeakage( 'dim9A.dat', 'leak9A.dat' )
     # simple 1-section polarizer

  a1 = 15.
  a2 = 60.
  dumpDimensions( 'dim9B.dat', r, a1, f, L180, a2, f, L90, 0., 0., 0.)
  computeLeakage( 'dim9B.dat', 'leak9B.dat' )
     # 2-section polarizer, maximally flat, from Kovac

  a1 = 15.
  a2 = 59.5
  dumpDimensions( 'dim9C.dat', r, a1, f, L180, a2, f, L90, 0., 0., 0.)
  computeLeakage( 'dim9C.dat', 'leak9C.dat' )
    # 2-section polarizer, CARMA design

  a1 = 6.05
  a2 = 28.63
  a3 = 67.59
  dumpDimensions( 'dim9D.dat', r, a1, f, L180, a2, f, L180, a3, f, L90 )
  computeLeakage( 'dim9D.dat', 'leak9D.dat' )
    # 3-section polarizer, maximally flat, from Kovac 

  a1 = 6.50
  a2 = 28.07
  a3 = 66.57
  dumpDimensions( 'dim9E.dat', r, a1, f, L180, a2, f, L180, a3, f, L90 )
  computeLeakage( 'dim9E.dat', 'leak9E.dat' )
    # wideband 3-section polarizer derived by Pancharatnam

def fig10() :
  r = .0235
  f = 0.006
  L90 = length( r=r, f=f, dphi0=90., fGHz=230. ) 
  L180 = 2. * L90
  a1 = 15.
  a2 = 60.
    # nominal dimensions of 2-section polarizer, maximally flat, from Kovac

  dumpDimensions( 'dim10a.dat', r, a1, f-.0002, L180, a2, f-.0002, L90, 0., 0., 0. )
  dumpDimensions( 'dim10a.dat', r, a1, f,       L180, a2, f,       L90, 0., 0., 0. )
  dumpDimensions( 'dim10a.dat', r, a1, f+.0002, L180, a2, f+.0002, L90, 0., 0., 0. )
  computeLeakage( 'dim10a.dat', 'leak10a.dat' )
    # change facet depths for both sections

  dumpDimensions( 'dim10b.dat', r, a1, f, L180, a2-1., f, L90, 0., 0., 0.)
  dumpDimensions( 'dim10b.dat', r, a1, f, L180, a2,    f, L90, 0., 0., 0.)
  dumpDimensions( 'dim10b.dat', r, a1, f, L180, a2+1., f, L90, 0., 0., 0.)
  computeLeakage( 'dim10b.dat', 'leak10b.dat' )
    # change angle between section 1 and section 2
  
# ---------------------------------------------------------------------------------------------------------- #
# return dimension = (xtarg +/- random error)
# error is gaussian distributed with sigma=xtol, but cannot exceed xtol
# ---------------------------------------------------------------------------------------------------------- #
def dimension( xtarg, xtol ) :
  if (xtarg == 0.) : return 0.
  x = xtarg + 2.*xtol                  # enter loop with dummy value guaranteed to be unacceptable
  while (abs(x - xtarg) > xtol) :
    x = random.gauss( xtarg, xtol )
  return x
   
# ---------------------------------------------------------------------------------------------------------- #
# writes file of polarizer dimensions with truncated Gaussian deviations (one line per polarizer) 
# ---------------------------------------------------------------------------------------------------------- #
def gaussdim( ntrials=200, dimfile='dims.txt', rtarg=0.0235, rtol=0.0001, ftarg=0.003, ftol=0.0001, 
              a1targ=15.0, a1tol=0.2, L1targ=0.4718, Ltol=0.001, a2targ=59.5, atol=0.1, L2targ=0.2359, 
              a3targ=0., L3targ=0. ) :
  ofile = open(dimfile, "w")
  
  # --- dump input parameters for future reference --- #
  ofile.write("# rtarg = %.6f, rtol = %.6f\n" % (rtarg,rtol) )
  ofile.write("# ftarg = %.6f, ftol = %.6f\n" % (ftarg,ftol) )
  ofile.write("# a1targ = %.2f, a1tol = %.2f, L1targ= %.6f, L1tol = %.6f \n" % (a1targ,a1tol,L1targ,Ltol) )
  ofile.write("# a2targ = %.2f, a2tol = %.2f, L2targ= %.6f, L2tol = %.6f \n" % (a2targ,atol,L2targ,Ltol) )
  ofile.write("# a3targ = %.2f, a3tol = %.2f, L3targ= %.6f, L3tol = %.6f \n" % (a3targ,atol,L3targ,Ltol) )

  # --- write one line per polarizer --- #
  for n in range(ntrials) :
    ofile.write(" %6.3f" % (1000.*dimension(rtarg,rtol)) )    # circ waveguide radius in mils
    ofile.write(" %7.2f" % dimension(a1targ,a1tol) )		  # angle of section 1 relative to input ref axis
    ofile.write(" %5.3f" % (1000.*dimension(ftarg,ftol)) )	  # facet depth of section 1, mils
    ofile.write(" %6.2f" % (1000.*dimension(L1targ,Ltol)) )   # length of section 1, mils
    ofile.write(" %7.2f" % dimension(a2targ,atol) )           # angle of section 2 relative to input section 1
    ofile.write(" %5.3f" % (1000.*dimension(ftarg,ftol)) )    # facet depth of section 2, mils  
    ofile.write(" %6.2f" % (1000.*dimension(L2targ,Ltol)) )   # length of section 2, mils
    ofile.write(" %7.2f" % dimension(a3targ,atol) )           # angle of section 3 relative to input section 2
    ofile.write(" %5.3f" % (1000.*dimension(ftarg,ftol)) )    # facet depth of section 3, mils  
    ofile.write(" %6.2f" % (1000.*dimension(L3targ,Ltol)) )   # length of section 3, mils
    ofile.write("\n")                             
  ofile.close()

# ---------------------------------------------------------------------------------------------------------- #
# leakages of polarizers with fabrication tolerances
# ---------------------------------------------------------------------------------------------------------- #
def fig11() :
  r = .0235
  f = 0.003
  L90 = length( r=r, f=f, dphi0=90., fGHz=230. ) 
  L180 = 2. * L90
  a1 = 15.
  a2 = 59.5
  gaussdim( ntrials=200, dimfile='dim11a.dat', rtarg=r, rtol=0.00010, ftarg=f, ftol=0.00010, 
      a1targ=a1, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=a2, atol=0.2, L2targ=L90 ) 
  gaussdim( ntrials=200, dimfile='dim11b.dat', rtarg=r, rtol=0.00015, ftarg=f, ftol=0.00015, 
      a1targ=a1, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=a2, atol=0.2, L2targ=L90 ) 
  gaussdim( ntrials=200, dimfile='dim11c.dat', rtarg=r, rtol=0.00020, ftarg=f, ftol=0.00020, 
      a1targ=a1, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=a2, atol=0.2, L2targ=L90 ) 
  computeLeakage( 'dim11a.dat', 'leak11a.dat' )
  computeLeakage( 'dim11b.dat', 'leak11b.dat' )
  computeLeakage( 'dim11c.dat', 'leak11c.dat' )

  f = 0.006
  L90 = length( r=r, f=f, dphi0=90., fGHz=230. ) 
  L180 = 2. * L90
  gaussdim( ntrials=200, dimfile='dim11d.dat', rtarg=r, rtol=0.00010, ftarg=f, ftol=0.00010, 
      a1targ=a1, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=a2, atol=0.2, L2targ=L90 ) 
  gaussdim( ntrials=200, dimfile='dim11e.dat', rtarg=r, rtol=0.00015, ftarg=f, ftol=0.00015, 
      a1targ=a1, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=a2, atol=0.2, L2targ=L90 ) 
  gaussdim( ntrials=200, dimfile='dim11f.dat', rtarg=r, rtol=0.00020, ftarg=f, ftol=0.00020, 
      a1targ=a1, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=a2, atol=0.2, L2targ=L90 ) 
  computeLeakage( 'dim11d.dat', 'leak11d.dat' )
  computeLeakage( 'dim11e.dat', 'leak11e.dat' )   # also used for fig 12a
  computeLeakage( 'dim11f.dat', 'leak11f.dat' )

def fig12b() :
  r = .0235
  f = 0.006
  L90 = length( r=r, f=f, dphi0=90., fGHz=230. ) 
  L180 = 2. * L90
  gaussdim( ntrials=200, dimfile='dim12b.dat', rtarg=r, rtol=0.00015, ftarg=f, ftol=0.00015, 
      a1targ=6.50, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=28.07, atol=0.2, L2targ=L180,
      a3targ=66.57, L3targ=L90) 
  computeLeakage( 'dim12b.dat', 'leak12b.dat' )

# ---------------------------------------------------------------------------------------------------------- #
# Fresnel formulae, transmission and reflection amplitude coeff, 1.5.2, eqn 20,21 (p. 40); for 1 surface
# ---------------------------------------------------------------------------------------------------------- #
def fresnel( n1, thetaI, n2, thetaT ) :
  tpar = 2. * n1 * math.cos(thetaI) / (n2 * math.cos(thetaI) + n1 * math.cos(thetaT))
  tperp = 2. * n1 * math.cos(thetaI) / (n1 * math.cos(thetaI) + n2 * math.cos(thetaT))
  rpar = (n2 * math.cos(thetaI) - n1 * math.cos(thetaT)) / (n2 * math.cos(thetaI) + n1 * math.cos(thetaT))
  rperp = (n1 * math.cos(thetaI) - n2 * math.cos(thetaT)) / (n1 * math.cos(thetaI) + n2 * math.cos(thetaT))
  return [tpar, tperp, rpar, rperp]

# ---------------------------------------------------------------------------------------------------------- #
# compute transmission through beamsplitter ("pellicle") of thickness tpel (inches)
# angI is the angle of incidence (degrees), nn is the refractive index of the beamsplitter material
# default nn=1.83 is appropriate for PETP = Mylar (Lamb 1996, Int. J. IR MM Waves, 17, pp. 1997-2034)
# returns  [tpar,tperp,Rpar,Rperp] 
# tpar and tperp are the AMPLITUDE TRANSMISSION coefficients (complex numbers)
# Rpar and Rperp are the POWER REFLECTION coefficients (real numbers)
# equations referenced to Born & Wolf 1970, Principles of Optics, Fourth Ed. (Oxford: Pergamon Press). 
# ---------------------------------------------------------------------------------------------------------- #
def pellicle( tpel=.001, angI=45, nn=1.83, fGHz=230. ) :
  if (tpel == 0.00) : return [1.,1.,0.,0.]
  thetaI = math.pi * angI / 180.	                        # convert to radians 
  thetaT = math.asin((math.sin(thetaI))/nn)	                        # Snell's law, eqn 8 (p. 38)
  [tpar,tperp,rpar,rperp] = fresnel( 1., thetaI, nn, thetaT )                        # entering dielectric
  [tprimepar,tprimeperp,rprimepar,rprimeperp] = fresnel( nn, thetaT, 1., thetaI )    # leaving dielectric
  delta = 4. * math.pi * 2.54 * tpel * nn * math.cos(thetaT) / (clight/fGHz)	  
    # phase shift of signal propagating once through the pellicle; eqn (1), p. 324
  result = []
  result.append(tpar*tprimepar/(1.-rpar*rpar*cmath.exp(1.j * delta)))
  result.append(tperp*tprimeperp/(1.-rperp*rperp*cmath.exp(1.j * delta)))
    # amplitude transmission coefficient, eqn (11), p. 325
  Fpar = 4.*rpar*rpar/pow((1. - rpar*rpar),2.)  
  Fperp = 4.*rperp*rperp/pow((1.-rperp*rperp),2.)
  sinsq = pow( math.sin(delta/2.), 2. )
  result.append( Fpar * sinsq/(1. + Fpar * sinsq) )	
  result.append( Fperp * sinsq/(1. + Fperp * sinsq) )
    # power reflection coefficients, eqn (15,16) p. 327 
  return result

# --- beamsplitter reflectivity --- #
def fig18( outfile='pelR.dat' ) :
  ofile = open( outfile, "w" )
  ofile.write("# tpel  Rpar35  Rperp35  Rpar45  Rperp45\n")
  for tpel in numpy.arange(0.00002,.00202,.00002) :
    [tpar,tperp,Rpar,Rperp] = pellicle( tpel=tpel, angI=35, fGHz=230.)
    ofile.write("%7.5f %10.4f %10.4f" % (tpel, 10.*math.log10(Rpar), 10.*math.log10(Rperp)))
    [tpar,tperp,Rpar,Rperp] = pellicle( tpel=tpel, angI=45, fGHz=230.)
    ofile.write("  %10.4f %10.4f\n" % (10.*math.log10(Rpar), 10.*math.log10(Rperp)))

# --- leakages including degradation by beamsplitter --- #
def fig19() :
  r = .0235
  f = 0.006
  L90 = length( r=r, f=f, dphi0=90., fGHz=230. ) 
  L180 = 2. * L90
  a1 = 15.
  a2 = 59.5
  dumpDimensions( 'nominal.dat', r, a1, f, L180, a2, f, L90, 0., 0., 0. )
  computeLeakage( 'nominal.dat', 'bs-none.dat', apel=0., tpel=0.000, aIpel=45. ) 
  computeLeakage( 'nominal.dat', 'bs-0-45-1mil.dat', apel=0., tpel=0.001, aIpel=45. ) 
  computeLeakage( 'nominal.dat', 'bs-45-45-1mil.dat', apel=45., tpel=0.001, aIpel=45. ) 
  computeLeakage( 'nominal.dat', 'bs-90-45-1mil.dat', apel=90., tpel=0.001, aIpel=45. ) 
  computeLeakage( 'dim11e.dat', 'leak11e-0-45-1mil.dat', apel=0., tpel=0.001, aIpel=45. ) 
  computeLeakage( 'dim11e.dat', 'leak11e-45-45-1mil.dat', apel=45., tpel=0.001, aIpel=45. ) 
  computeLeakage( 'dim11e.dat', 'leak11e-90-45-1mil.dat', apel=90., tpel=0.001, aIpel=45. ) 

# ---------------------------------------------------------------------------------------------------------- #
# compute leakage from HFSS simulation results
# normally the HFSS calculation launches an X-polarized signal, returns S21(X,Y) 
# each infile line: f_GHz, S21(X:X)amp_dB, S21(X:X)phs_degr, S21(X:Y)amp_dB, S21(X:Y)phs_degr
# ---------------------------------------------------------------------------------------------------------- #
def HFSSleak( infile, leakfile ) :
  infile = open( infile, "r" )
  ofile = open( leakfile, "w")
  for line in infile:
    if line.startswith("#"):
      ofile.write( line )
    else:
      a = line.split()                                # split line into string tokens
      fGHz =  float(a[0]) 	                      # frequency
      ampX = math.sqrt(pow(10.,float(a[1])/10.))      # convert amp to VOLTAGE
      phsX = math.pi * float(a[2]) / 180.	      # convert phs to RADIANS
      ampY = math.sqrt(pow(10.,float(a[3])/10.))
      phsY = math.pi * float(a[4]) / 180.
      vec = numpy.array( [ (ampX*math.cos(phsX) + 1j*ampX*math.sin(phsX)),
                     (ampY*math.cos(phsY) + 1j*ampY*math.sin(phsY)) ] )
      print vec, numpy.dot(vec,R), numpy.dot(vec,L)
      Rleak = abs(numpy.dot(vec, R)) 		# amplitude of single complex number 
      Lleak = abs(numpy.dot(vec, L))
      ofile.write("%8.3f %8.3f %8.3f %8.3f %8.3f %10.5f %10.5f\n" %
        (fGHz, ampX, float(a[2]), ampY, float(a[4]), Rleak, Lleak))
  ofile.close()
  
# ---------------------------------------------------------------------------------------------------------- #
# analytic simulation of full polarizer including curved transitions
# ---------------------------------------------------------------------------------------------------------- #
def fig16( outfile='fullsim1.dat' ) :
  ofile = open( outfile, "w")
  for freq in numpy.arange(200.,271.,1.) :
    v1 = Y
    v2 = Jrot( v1, -15. )
    [dx, ph] = dphitaper2 ( r=0.0235, f=0.006, Rc=0.125, fGHz=freq, nsteps=1000 ) 
    phshift = 2.*ph + dphi( r=0.0235, f=0.006, L=0.1058, fGHz=freq ) 
    v3 = Jdelay( v2, phshift )
    v2 = Jrot( v3, -59.5 )
    phshift = 2.*ph + dphi( r=0.0235, f=0.006, L=0.0302, fGHz=freq ) 
    v3 = Jdelay( v2, phshift )
    v2 = Jrot( v3, 74.5 )
    vout = v2
    Rleak = abs(numpy.dot(vout, R)) 		# amplitude of single complex number 
    Lleak = abs(numpy.dot(vout, L))
    ofile.write("%6.1f %8.5f %8.5f\n" % (freq, Rleak, Lleak ) )
  ofile.close()

# ==== superfluous older routines below this line ===== #

def finalcheck3( outfile='reversed.dat' ) :
  ofile = open( outfile, "w")
  for freq in numpy.arange(200.,271.,1.) :
    v1 = L
    v2 = Jrot( v1, 74.5 )
    phshift = dphi( r=0.0235, f=0.006, L=0.0752, fGHz=freq ) 
    v3 = Jdelay( v2, phshift )
    v2 = Jrot( v3, -59.5 )
    phshift = dphi( r=0.0235, f=0.006, L=0.1504, fGHz=freq ) 
    v3 = Jdelay( v2, phshift )
    vout = Jrot( v3, -15. )
    cleak = numpy.dot(vout, Y) 		# single complex number 
    ofile.write("%6.1f %6.4f\n" % (freq, abs(cleak)) )
  ofile.close()



# ---------------------------------------------------------------------------------------------------------- #
# for each set of (2-section) polarizer dimensions in infile, compute X- and Y- output powers in Perficio test
# in this test, linear polarized signal is injected at angle a3 into the circ pol output, and we
#   measure the relative powers from the X- and Y- outputs of an OMT that is connected to the other port
# note that, as usual, a1 is the angle of the halfwave section relative to the OMT reference axis, etc.
# ---------------------------------------------------------------------------------------------------------- #
def perflist( dimfile, outfile ) :
  r = []	# empty list of circular waveguide radii
  a1 = []   # angle of section 1 relative to OMT reference axis
  f1 = []   # facet depth of section 1
  L1 = []   # length of section 1
  a2 = []   # angle of section 2 relative to section 1
  f2 = []   # facet depth of section2
  L2 = []   # length of section 2
  a3 = []   # angle of input linear pol relative to section 2
  
  # --- read the dimensions into internal lists --- #
  infile = open(dimfile, "r")
  ofile = open(outfile, "w")
  ntrials = 0
  for line in infile:
    if line.startswith("#"):
      ofile.write(line)
    else:
      a = line.split()                            # split line into string tokens
      r.append( float(a[0])/1000. )		# convert from mils to inches
      a1.append( float(a[1]) )
      f1.append( float(a[2])/1000. )	    # convert facet depth from mils to inches
      L1.append( float(a[3])/1000. )
      a2.append( float(a[4]) )
      f2.append( float(a[5])/1000. )	    # convert facet depth from mils to inches
      L2.append( float(a[6])/1000. )
      a3.append( float(a[7]) )
      ntrials = ntrials + 1
  infile.close()
    
  # --- compute power from OMT X and Y ports for each polarizer in the list --- #
  freq = numpy.arange(205.,271.)
  powX = numpy.zeros([len(freq),ntrials],Float)
  powY = numpy.zeros([len(freq),ntrials],Float)
  for n in range(ntrials) :
    m = 0									# freq index
    for fGHz in freq :
      v1 = Y
      v2 = Jrot( v1, a3[n] )
      v3 = Jdelay( v2, dphi( r[n], f2[n], L2[n], fGHz ) )
      v2 = Jrot( v3, a2[n] )
      v3 = Jdelay( v2, dphi( r[n], f1[n], L1[n], fGHz ) )
      v2 = Jrot( v3, a1[n] )
      powX[m,n] = pow(abs(numpy.dot(v2,X)),2.)
      powY[m,n] = pow(abs(numpy.dot(v2,Y)),2.)
      m = m + 1

  # --- write dimensions to output file --- #
  dumplist ( ofile, r,  "# r  :", "    %6.3f   ", 1000.) 
  dumplist ( ofile, a1, "# a1 :", "    %6.2f   ", 1.) 
  dumplist ( ofile, f1, "# f1 :", "    %6.3f   ", 1000.) 
  dumplist ( ofile, L1, "# L1 :", "    %6.2f   ", 1000.) 
  dumplist ( ofile, a2, "# a2 :", "    %6.2f   ", 1.) 
  dumplist ( ofile, f2, "# f2 :", "    %6.3f   ", 1000.) 
  dumplist ( ofile, L2, "# L2 :", "    %6.2f   ", 1000.) 
  dumplist ( ofile, a3, "# a3 :", "    %6.2f   ", 1.) 

  # --- dump powers at each freq, for all polarizers; also compute avg and rms --- #
  m = 0
  for fGHz in freq :
    ofile.write("%6.1f" % fGHz)
    sum = 0.
    for n in range(ntrials) :
      ofile.write(" %6.4f %6.4f" % (powX[m,n], powY[m,n]) )
      sum = sum + powX[m,n]
    avg = sum/ntrials
    var = 0.
    for n in range(ntrials) :
      dif = powX[m,n] - avg
      var = var + dif*dif
    if (ntrials > 1) :
      var = var/(ntrials-1)
    ofile.write("  %6.4f %6.4f\n" % (avg,math.sqrt(var)))
    m = m + 1
  ofile.close()     


# --- simulate test data for Perficio polarizers --- #
# --- inject linear pol at angle1, compute emergent X and Y power --- #
def perf( angle1=15.5, angle2=59.5, angle3=15., fdepth=.006, outfile='perf.dat' ) :
  ofile = open( outfile, "w")
  ofile.write("# angle1 = %.2f  # incoming linear to quarter wave section\n" % angle1)
  ofile.write("# angle2 = %.2f  # quarter wave to halfwave\n" % angle2)
  ofile.write("# angle3 = %.2f  # halfwave to OMT\n" % angle3)
  ofile.write("# fdepth = %.4f  # facet depth\n" % fdepth)
  for freq in numpy.arange(205.,271.,1.) :
    v1 = Y
    v2 = Jrot( v1, angle1 )
    [dx, ph] = dphitaper2 ( r=0.0235, f=fdepth, Rc=0.0625, fGHz=freq, nsteps=100 ) 
    phshift1 = 2.*ph + dphi( r=0.0235, f=fdepth, L=0.0444, fGHz=freq ) 
    v3 = Jdelay( v2, phshift1 )
    v2 = Jrot( v3, angle2 )
    phshift2 = 2.*ph + dphi( r=0.0235, f=fdepth, L=0.1196, fGHz=freq ) 
    v3 = Jdelay( v2, phshift2 )
    v2 = Jrot( v3, angle3 )
    magX = abs(numpy.dot(v2,X))
    magY = abs(numpy.dot(v2,Y))
    # --- also compute leakage --- #
    v1 = R
    v2 = Jrot( v1, angle1 )
    v3 = Jdelay( v2, phshift1 )
    v2 = Jrot( v3, angle2 )
    v3 = Jdelay( v2, phshift2 )
    v2 = Jrot( v3, angle3 )
    cleak = numpy.dot(v2, X) 		# single complex number 
   
    ofile.write("%6.1f %6.4f %6.4f %6.4f\n" % (freq, pow(magX,2.), pow(magY,2.), abs(cleak)) )
  ofile.close()


def fig2 ( ) :
  ofile = open("fig2.dat","w")
  for f1 in numpy.arange(0.0020, 0.0071, .0001):
    L1 = length( f=f1 )
    ofile.write("%6.4f %6.4f %6.2f\n" % ( f1, L1, dphi( f=(f1+.0001), L=L1 ) ))
  ofile.close()

# ---------------------------------------------------------------------------------------------------------- #
# writes file of polarizer dimensions (one line per polarizer), stepping one or more over linear range 
# rtarg = radius of circular guide; ftarg = facet depth; a1targ,a2targ,a3targ = angles of sections 1,2,3;
# L1targ,L2targ,L3targ = lengths of sections 1,2,3
# ---------------------------------------------------------------------------------------------------------- #
def lindim( ntrials=1, dimfile='dims.txt', rtarg=0.0235, rstep=0.0, ftarg=0.006, fstep=0.0, 
            a1targ=15.0, a1step=0.0, L1targ=0.1658, L1step=0.0,
            a2targ=59.5, a2step=0.0, L2targ=0.0829, L2step=0.0,
            a3targ=0., a3step=0.0, L3targ=0., L3step=0.0 ) :
  ofile = open(dimfile, "w")
  ofile.write("# rtarg = %.6f, rstep = %.6f\n" % (rtarg,rstep) )
  ofile.write("# ftarg = %.6f, fstep = %.6f\n" % (ftarg,fstep) )
  ofile.write("# a1targ = %.2f, a1step = %.2f, L1targ= %.6f, L1step = %.6f \n" % (a1targ,a1step,L1targ,L1step) )
  ofile.write("# a2targ = %.2f, a2step = %.2f, L2targ= %.6f, L2step = %.6f \n" % (a2targ,a2step,L2targ,L2step) )
  ofile.write("# a3targ = %.2f, a3step = %.2f, L3targ= %.6f, L3step = %.6f \n" % (a3targ,a3step,L3targ,L3step) )
  ofile.close()
  nn = ntrials/2

  # --- write out nominal design first - guarantees that there is at least one output line --- #
  dumpDimensions( dimfile, rtarg, a1targ, ftarg, L1targ, a2targ, ftarg, L2targ, a3targ, ftarg, L3targ, a4targ )

  # --- now write out one line per step for variables with step .ne. 0 --- #
  if (rstep != 0.) :
    for r in numpy.arange( rtarg-nn*rstep, rtarg+(nn+1)*rstep, rstep ) :
      dumpDimensions( dimfile, r, a1targ, ftarg, L1targ, a2targ, ftarg, L2targ, a3targ, ftarg, L3targ )
  if (fstep != 0.) :
    for f in numpy.arange( ftarg-nn*fstep, ftarg+(nn+1)*fstep, fstep ) :
      dumpDimensions( dimfile, rtarg, a1targ, f, L1targ, a2targ, f, L2targ, a3targ, f, L3targ, a4targ )
  if (a1step != 0.) :
    for a1 in numpy.arange( a1targ-nn*a1step, a1targ+(nn+1)*a1step, a1step ) :
      dumpDimensions( dimfile, rtarg, a1, ftarg, L1targ, a2targ, ftarg, L2targ, a3targ, ftarg, L3targ )
  if (L1step != 0.) :
    for L1 in numpy.arange( L1targ-nn*L1step, L1targ+(nn+1)*L1step, L1step ) :
      dumpDimensions( dimfile, rtarg, a1targ, ftarg, L1, a2targ, ftarg, L2targ, a3targ, ftarg, L3targ )
  if (a2step != 0.) :
    for a2 in numpy.arange( a2targ-nn*a2step, a2targ+(nn+1)*a2step, a2step ) :
      dumpDimensions( dimfile, rtarg, a1targ, ftarg, L1targ, a2, ftarg, L2targ, a3targ, ftarg, L3targ )
  if (L2step != 0.) :
    for L2 in numpy.arange( L2targ-nn*L2step, L2targ+(nn+1)*L2step, L2step ) :
      dumpDimensions( dimfile, rtarg, a1targ, ftarg, L1targ, a2targ, ftarg, L2, a3targ, ftarg, L3targ )
  if (a3step != 0.) :
    for a3 in numpy.arange( a3targ-nn*a3step, a3targ+(nn+1)*a3step, a3step ) :
      dumpDimensions( dimfile, rtarg, a1targ, ftarg, L1targ, a2targ, ftarg, L2targ, a3, ftarg, L3targ )
  if (L3step != 0.) :
    for L3 in numpy.arange( L3targ-nn*L3step, L3targ+(nn+1)*L3step, L3step ) :
      dumpDimensions( dimfile, rtarg, a1targ, ftarg, L1targ, a2targ, ftarg, L2targ, a3targ, ftarg, L3 )

# ---------------------------------------------------------
# returns cutoff freqs for circular waveguide of radius r inches
#  from Harvey pp 16-17: lamdaCutoff/r = 2.6127 for TM01, 
#    1.6398 for TM11, TE01, 3.4126 for TE11, 2.0572 for TE21
# ----------------------------------------------------------
def modeCutoffs ( r=.0235 ) :
  print "diam = %.4f inches:" % (2.*r)
  r = 2.54 * r
  c = 29.9792458
  print "  %.4f GHz TE11" % (c/(3.4126 * r))
  print "  %.4f GHz TM01" % (c/(2.6127 * r))
  print "  %.4f GHz TE21" % (c/(2.0572 * r))
  print "  %.4f GHz TM11,TE01" % (clight/(1.6398 * r))
  fratio = 220./(c/(3.4126*r)) # f/fc
  alphac = (0.423/pow(r,1.5)) * \
    ((pow(fratio,-0.5) + pow(fratio,1.5)/2.38)/(math.sqrt(pow(fratio,2.) - 1.)))
  print "  %.2f dB/inch loss at 220 GHz" % (alphac/1200.)

 
# --- compute tsys expected for beamsplitter transmission trans ---
def tsys( trcvr, trans ) : 
  tcold = 77.
  tamb = 297.
  tsys = (tamb-tcold)*(trcvr + trans*tcold + (1.-trans)*tamb)/(tamb - trans*tcold - (1.-trans)*tamb) - tcold
  return tsys 

def perfplots() :
  r = .0235
  f = 0.006
  a1 = 15.
  a2 = 59.5
  L90 = length( r=r, f=f, dphi0=90., fGHz=230. ) 
  L180 = 2. * L90

  #gaussdim( ntrials=200, dimfile='dim1.dat', rtarg=r, rtol=0.00015, ftarg=f, ftol=0.00015, 
  #    a1targ=a1, a1tol=0.2, L1targ=L180, Ltol=0.001, a2targ=a2, atol=0.2, L2targ=L90, a3targ=15. ) 
  #perflist("dim1.dat","meetsSpec.dat")		

  #gaussdim( ntrials=200, dimfile='dim2.dat', rtarg=r, rtol=0.0005, ftarg=f, ftol=0.0008, 
  #    a1targ=a1, a1tol=2, L1targ=L180, Ltol=0.005, a2targ=a2, atol=2, L2targ=L90, a3targ=15. ) 
  #perflist("dim2.dat","badSpec.dat")		

  # gaussdim( ntrials=1, dimfile='dim3.dat', rtarg=r, rtol=0.000, ftarg=f, ftol=0.000, 
  #    a1targ=a1, a1tol=0, L1targ=L180, Ltol=0.000, a2targ=a2, atol=0, L2targ=L90, a3targ=15. ) 
  #perflist("dim3.dat","perfect.dat")		

  # --- pretty good match to the test results --- #
  #  gaussdim( ntrials=200, dimfile='dim4.dat', rtarg=r, rtol=0.0005, ftarg=f, ftol=0.0008, 
    #    a1targ=a1, a1tol=2., L1targ=L180, Ltol=0.005, a2targ=a2, atol=2, L2targ=L90, a3targ=15. )
    #computeLeakage( dimfile='dim4.dat', leakfile='leak4.dat', efficfile='effic4.dat' ) 
    #perflist("dim4.dat","dim4Spec.dat")		

  gaussdim( ntrials=200, dimfile='dim5.dat', rtarg=r, rtol=0.00015, ftarg=f, ftol=0.00015, 
      a1targ=a1, a1tol=0.2, L1targ=L180-.04, Ltol=0.001, a2targ=70., atol=0.2, L2targ=L90, a3targ=15. )
  computeLeakage( dimfile='dim5.dat', leakfile='leak5.dat', efficfile='effic5.dat' ) 


# --- emulate tests of CP11 done at Agilent on 29jan2010 --- #
def Agilentcheck( outfile='Agilent1.dat' ) :
  ofile = open( outfile, "w")
  for freq in numpy.arange(200.,271.,1.) :
    v1 = Y
    v2 = Jrot( v1, -15. )
    [dx, ph] = dphitaper2 ( r=0.0235, f=0.006, Rc=0.0625, fGHz=freq, nsteps=100 ) 
    phshift = 2.*ph + dphi( r=0.0235, f=0.006, L=0.1196, fGHz=freq ) 
    print freq, ph, phshift
    v3 = Jdelay( v2, phshift )
    v2 = Jrot( v3, -59.5 )
    phshift = 2.*ph + dphi( r=0.0235, f=0.006, L=0.0444, fGHz=freq ) 
    v3 = Jdelay( v2, phshift )
    v2 = Jrot( v3, 74.5 )
    Xout = numpy.dot(v2, X) 		# single complex number 
    Yout = numpy.dot(v2, Y)
    phsX = 180. * math.atan2(Xout.imag,Xout.real) / math.pi
    phsY = 180. * math.atan2(Yout.imag,Yout.real) / math.pi
    ofile.write("%6.1f %6.4f %8.4f %6.4f %8.4f %10.4f\n" % (freq, abs(Xout), phsX, abs(Yout), phsY, phsX-phsY)) 
  ofile.close()

# ---------------------------------------------------------------------------------------------------------- #
# compute differential phase shift through linear tapered section at frequency fGHz
# the facet depth starts at 0.000 and increases linearly to f; all dimensions in inches
# ---------------------------------------------------------------------------------------------------------- #
def dphitaper( r=0.0235, f=0.006, L=.04, fGHz=230., nsteps=1000 ) :
  dL = L/nsteps      # increment in x
  df = f/nsteps      # increment in facet depth
  x = dL/2.          # mean x-coordinate of 1st section
  fx = f - df/2.     # mean facet depth of 1st section
  totphase = 0.
  for n in range(nsteps) :
    totphase = totphase + dphi(r=r, f=fx, L=dL, fGHz=fGHz)
    fx = fx - df
  return totphase

# ---------------------------------------------------------------------------------------------------------- #
# compute optimum facet depth and length for stepped matching section from circular to faceted guide
# ---------------------------------------------------------------------------------------------------------- #
def stepmatch(dfacet=.006, r=0.0235, fGHz=230., outfile="match.dat" ) :
  ofile = open(outfile,"w")
  [fcx0,fcy0] = cutoff( r, 0. )
  [fcx1,fcy1] = cutoff( r, dfacet/r )
  Zx0 = fGHz * 377. / math.sqrt(fGHz*fGHz - fcx0*fcx0)
  Zy0 = fGHz * 377. / math.sqrt(fGHz*fGHz - fcy0*fcy0)
  Zx1 = fGHz * 377. / math.sqrt(fGHz*fGHz - fcx1*fcx1)
  Zy1 = fGHz * 377. / math.sqrt(fGHz*fGHz - fcy1*fcy1)
  vRx = (Zx1 - Zx0) / (Zx1 + Zx0)
  vRy = (Zy1 - Zy0) / (Zy1 + Zy0)
  print "voltage, power reflection coeff with no matching section"
  print " .. x: %.3e %.2f" % (vRx, -20.*log(abs(vRx))/2.30259 )
  print " .. y: %.3e %.2f" % (vRy, -20.*log(abs(vRy))/2.30259 )
  ofile.write("# dmatch   fcutoffX    Zx-Zxm    fcutoffY    Zy-Zym\n") 
  Zxm = math.sqrt(Zx0*Zx1)
  Zym = math.sqrt(Zy0*Zy1)
  dmatchx = 0.
  dmatchy = 0.
  Zdifxleast = 100000.
  Zdifyleast = 100000.
  for d in numpy.arange (0., dfacet, .00001) :
    [fcxm,fcym] = cutoff( r, d/r )
    Zx = fGHz * 377. / math.sqrt(fGHz*fGHz - fcxm*fcxm)
    Zy = fGHz * 377. / math.sqrt(fGHz*fGHz - fcym*fcym)
    ofile.write("%8.5f  %10.3f %10.6f  %10.3f %10.6f\n" % (d, fcxm,(Zx-Zxm),fcym,(Zy-Zym))) 
    if (abs(Zx-Zxm) < Zdifxleast) :
       dmatchx = d
       Zdifxleast = abs(Zx-Zxm)
    if (abs(Zy-Zym) < Zdifyleast) :
       dmatchy = d
       Zdifyleast = abs(Zy-Zym)
  [fcxm,fcym] = cutoff( r, dmatchx/r )
  lmatch = clight/(2.54 * 4. * math.sqrt(fGHz*fGHz - fcxm*fcxm)) 		# lambda_guide/4.
  print "best X match is depth %.5f, length %.4f" % (dmatchx,lmatch)
  [fcxm,fcym] = cutoff( r, dmatchy/r )
  lmatch = clight/(2.54 * 4. * math.sqrt(fGHz*fGHz - fcym*fcym)) 
  print "best Y match is depth %.5f, length %.4f" % (dmatchy,lmatch)
  ofile.close()

# --- check fresnel by reproducing Fig 1.12 in Born and Wolf --- #
def BWfig12():
  ofile = open( "BWfig12.dat", "w" )
  nn = 1.52
  for angI in numpy.arange(0., 90., .5) :
    thetaI = math.pi * angI / 180.	  
    thetaT = math.asin((math.sin(thetaI))/nn)
    [tpar,tperp,rpar,rperp] = fresnel( 1., thetaI, nn, thetaT ) 
    ofile.write("%.2f %.4f %.4f\n" % (angI, rpar*rpar, rperp*rperp))
  ofile.close()
    
def testpellicle(ang) :
  ofile = open("testpellicle", "w")
  for tp in numpy.arange(0.,.05,.0001) :
      [tpar,tperp,rpar,rperp] = pellicle( tpel=tp, angI=ang )
      ofile.write("%6.4f %8.5f %8.5f\n" % (tp,rpar,rperp))
  ofile.close()

# --- added 18apr2010 - check whether DR and DL can have different magnitudes --- #
def LeakageRealityCheck() :
  v1 = R
  v2 = Jrot(v1, -15.5)
  v3 = Jdelay(v2, 85)
  v2 = Jrot(v3, -60.)
  v3 = Jdelay(v2, 180)
  v2 = Jrot(v3, -15.)
  print "D_L ", numpy.dot(v2,X), abs(numpy.dot(v2,X))   # note: R* = L
  print "D_R ", numpy.dot(v2,Y), abs(numpy.dot(v2,Y))   # note: L* = R
  v1 = L
  v2 = Jrot(v1, -15.5)
  v3 = Jdelay(v2, 85)
  v2 = Jrot(v3, -60.)
  v3 = Jdelay(v2, 180)
  v2 = Jrot(v3, -15.)
  print "D_L ", numpy.dot(v2,X), abs(numpy.dot(v2,X))   # note: R* = L
  print "D_R ", numpy.dot(v2,Y), abs(numpy.dot(v2,Y))   # note: L* = R

# --- added 21jan2015 - test whether DL = DR* for simple pol
def LeakageRealityCheck2() :
  v1 = R
  v2 = Jrot(v1, -45)
  v3 = Jdelay(v2, 85)
  v2 = Jrot(v3, 45.)
  print "D_L ", numpy.dot(v2,X), abs(numpy.dot(v2,X))   # note: R* = L
  print "D_R ", numpy.dot(v2,Y), abs(numpy.dot(v2,Y))   # note: L* = R
  v1 = L
  v2 = Jrot(v1, -45)
  v3 = Jdelay(v2, 85)
  v2 = Jrot(v3, 45.)
  print "D_L ", numpy.dot(v2,X), abs(numpy.dot(v2,X))   # note: R* = L
  print "D_R ", numpy.dot(v2,Y), abs(numpy.dot(v2,Y))   # note: L* = R


def plotVecs( infile ) :
  DLx = []
  DLy = []
  DRx = [] 
  DRy = []
  max = 0.1
  #py.ion()
  #py.c4lf()
  py.axis( [-max,max,-max,max] )
  fin = open( infile, "r" )
  for line in fin:
    if not line.startswith("#") :
      a = line.split()
      DLx.append(  a[2] )
      DLy.append(  a[3] )
      DRx.append(  a[5] )
      DRy.append(  a[6] )
      py.plot( numpy.array(DLx), numpy.array(DLy), 'rD', markersize=6., linewidth=4 )
      py.plot( numpy.array(DRx), numpy.array(DRy), 'b+', markersize=6., linewidth=4 )
  fin.close() 
  py.xlabel('real')
  py.ylabel('imag')
  py.grid(True)
  py.axes().set_aspect('equal')
  # py.figlegend( ('r+','g+'), ('Dx','Dy'), 'upper right')   saw nothing
  py.text( -0.95*max, .9*max, infile )
  #py.annotate( filename, (0., 0.) )
  #py.savefig("testfig.ps")    
  py.show()
  #py.draw()
    

# ---------------------------------------------------------------------------------------------------------- #
# propagate Y, R, L from antenna to OMT through beamsplitter and 2-section polarizer to simulate xyauto and 
#    leakage measurements
# assume r = 0.0235 circular waveguide radius
# X-axis is along the meridian ('vert' pol), Y-axis is parallel to the horizon ('horiz' pol)
#    
#    tpel = beamsplitter ("pellicle") thickness in inches
#    aIpel = angle of incidence to the beamsplitter, degrees (30 for bima, 45 for ovro)
#    aPpel = angle from X-axis to plane of incidence of beamsplitter, degrees (90 for bima, ? for ovro)
#    f = facet depth, inches (nominal = 0.006)
#    L180 = length of halfwave section (nominal = 0.15126)
#    L90 = length of quarterwave section (nominal = 0.07563)
#    a1 = angle of halfwave section relative to OMT Y-axis (nominal = 15)
#    a2 = angle of quarterwave section relative to halfwave section (nominal = 59.5)
#
# comment added 23jun2014: applying xyauto as done in this routine makes no sense because the real
#    xyauto folds together the 2 sidebands and takes the average
# ----------------------------------------------------------------------------------------------------------- #
def computeLeakage2( tpel=0., aIpel=30., aPpel=90., fdepth=0.006, a1=15., a2=59.5, L90=0.07563, \
     L180=0.15126, outfile="junk", xyphase=False) :
  r = 0.0235
  freqsave = []
  phdifsave = []
  ofile = open(outfile, "w")
  ofile.write("# beamsplitter thickness tpel = %.4f\n" % tpel)
  ofile.write("# beamsplitter angIncidence aIpel = %.2f\n" % aIpel)
  ofile.write("# beamsplitter angPlaneIncidence aPpel = %.2f\n" % aPpel)
  ofile.write("# facet depth fdepth = %.4f\n" % fdepth)
  ofile.write("# L90 = %.5f\n" % L90)
  ofile.write("# L180 = %.5f\n" % L180)
  ofile.write("# angle of halfwave section relative to OMT X-axis a1 = %.2f\n" % a1)
  ofile.write("# angle of quarterwave section relative to halfwave section a2 = %.2f\n" % a2)
  ofile.write("#  fGHz  R_mag  R_real  R_imag    L_mag  L_real L_imag   XampA  XphsA  XampB  XphsB   phsdifAB\n" )
  veclist = [X,R,L]
  for fGHz in numpy.arange(200.,271.,1.) :
    print " "
    ofile.write("%6.1f" % fGHz)
    phshift1 = dphi( r=0.0235, f=fdepth, L=L180, fGHz=fGHz )
    phshift2 = dphi( r=0.0235, f=fdepth, L=L90, fGHz=fGHz )
    nvec = 0
    for vec in veclist :

      v1 = Jrot( vec, aPpel)                   # rotate into plane of incidence of beamsplitter
      v2 = Jbsplit( v1, tpel, aIpel, fGHz )	   # propagate X,Y through beamsplitter
      v1 = Jrot( v2, -1.*aPpel )	           # rotate back to original reference frame
      v2 = Jrot( v1, -1.*(a1+a2) )             # angle of quarterwave section relative to telescope X-axis
      v1 = Jdelay( v2, phshift2 )              #  delay through quarterwave section
      v2 = Jrot( v1, a2 )			           # angle of halfwave relative to quarterwave section
      v1 = Jdelay( v2, phshift1 )	           # delay through halfwave section
      v2 = Jrot( v1, a1 )			           # rotate into OMT frame

      [ampA,phsA,ampB,phsB,phsdifAB] = ampPhs(v2)
      print "%5.0f   %7.5f %7.2f %7.5f %7.2f  %7.2f" % (fGHz, ampA,phsA,ampB,phsB,phsdifAB)
        # ... print amp,phs of X and Y components, and phs(X)-phs(Y)

      if (nvec == 0) :				
        [ampX,phsX,ampY,phsY,phsdif] = ampPhs(v2)	 # save values for noise source to print at end of line
        freqsave.append(fGHz)
        phdifsave.append(phsdif)
      else :
        leak = numpy.dot(v1,X)
        if (abs(leak) > abs(numpy.dot(v1,Y))) :
          leak = numpy.dot(v1,Y)
        ofile.write("     %7.5f %8.5f %8.5f" % (abs(leak),leak.real,leak.imag))
        if (nvec == 2 ) :							# print noise source values at end of line
          ofile.write("    %7.5f %7.2f %7.5f %7.2f  %7.2f" % (ampX,phsX,ampY,phsY,phsdif))
      nvec = nvec + 1
    ofile.write("\n")

  ofile.close()

  for n in range( 8, len(freqsave)-8 ) :
    pLSB = phdifsave[n] - phdifsave[n-8]
    if pLSB > 300. :
      pLSB = pLSB - 360.
    if pLSB < -300. :
      pLSB = pLSB + 360.
    pUSB = phdifsave[n+8] - phdifsave[n]
    if pUSB > 300. :
      pUSB = pUSB - 360.
    if pUSB < -300. :
      pUSB = pUSB + 360.
    print freqsave[n-8],freqsave[n],freqsave[n+8],pLSB, pUSB, pLSB-pUSB

  #plotVecs( outfile )


def doit() :
  #computeLeakage3( angle1=-45., angle3=45., outfile="aleak.1") 
  #computeLeakage3( angle1=-47., angle3=47., outfile="aleak.2") 
  #computeLeakage3( angle1=-47., angle3=45., outfile="aleak.3") 
  tpel = 0.
  aIpel = 30.
  #computeLeakage2( tpel, aIpel, fdepth=0.006, angle1=-74.5, angle2=59.5, angle3=15.0, outfile="aleak.4") 
  #computeLeakage2( 0., 30., fdepth=0.006, angle1=0., angle2=59.5, angle3=14.0, outfile="aleak.5") 
  # below: 22jun2014, try to figure out x-y phase
  print "begin computeLeakage2"
  computeLeakage2( tpel, aIpel, outfile="xyautoSimulate.dat", xyphase=False) 

# ============== add routines from waveplate.py below ================= #


# 5plates: 7, 36, 102, 36, 7 degrees
  
def nsapphire( fGHz ) :
  # equations from Savini et al. 2006, Appl Optics, p 8912 do not recover numbers
  #   that they give in Fig 4
  # no = 3.053 + (4.7e-4 * fGHz) + (2.2e-10 * pow(fGHz,2.)) + (1.1e-12 * pow(fGHz,3.))
  # ne = 3.387 + 1.3e-5 * fGHz
  # instead, use values from Afsar 1987, IEEE-MTT, Fig. 3
  no = 3.0647
  ne = 3.4037
  return [no,ne]

# effective index for sapphire with polarization at angle theta to ordinary axis
# e.g., if polarization is perpendicular to theta, neff = no; if parallel, neff = ne
def neff( thetaDeg ) :
  theta = math.pi*thetaDeg/180.
  no,ne = nsapphire( 150. )    # take 150 GHz as center freq; it doesn't matter!
  nsq = pow( no*ne, 2. ) / ( pow(no*math.cos(theta),2.) + pow(ne*math.sin(theta),2.) )
  return math.sqrt(nsq)

def doit6() :
  stack2 = [ [22.5, .37], [22.5, .37] ]
  stack3 = [ [15., .37], [75., .37], [15.,.37]]
  #v2 = propagate( 150., X, stack1 )
  makepath( 150., X, stack1 )
  #print "Stokes for X: ", stokes(X)
  #print "Stokes for Y: ", stokes(Y)
  #print "Stokes for R: ", stokes(R)
  #print "Stokes for L: ", stokes(L)
  #print "Stokes for rotated R: ", stokes( Jrot( R, 45 ))

# propagate polarization state v1 through waveplate 'plate'
# v1 is a complex matrix in the X and Y basis
# plate is a list of angles (relative to the X axis) and thicknesses
def propagate( fGHz, vin, stack ) :
  v1 = vin
  for plate in stack :
    print "propagate through plate with angle %.2f degrees, thickness %.3f cm" % (plate[0],plate[1])
    v2 = Jrot( v1, plate[0] )
    v3 = Jdelay( v2, dphi( plate[1], fGHz ) )
    v1 = Jrot( v3, -plate[0])
  return v1
  

# Stokes Q,U,V are x,y,z of the poincare sphere!
def stokes( vin ) :
  Ex = vin[0]
  Ey = vin[1]
  Q = (Ex * Ex.conjugate() - Ey * Ey.conjugate())
  U = (Ex * Ey.conjugate() + Ex.conjugate() * Ey)
  V = (cmath.exp(-1j * math.pi/2.) * ( Ex * Ey.conjugate() - Ex.conjugate() * Ey ))
  return numpy.real( [Q,U,V] )

# poincare plot

# generates array [x,y,z] = [Q,U,V] showing path from start to final 
def makepath( fGHz, vin, stack ) :
  path = []
  v1 = vin
  for plate in stack :
    print v1
    v2 = Jrot( v1, plate[0] )                   # v2 is rotated starting point for this segment
    for t in numpy.arange( 0., plate[1], .005 ) :
      v3 = Jdelay( v2, dphi(L=t, fGHz=fGHz) )      # for sapphire put in different dphi!
      v4 = Jrot( v3, -plate[0] )	            # v4 is endpoint for thickness t
      path.append( stokes( v4 ) )
    v3 = Jdelay( v2, dphi( L=plate[1], fGHz=fGHz ) )  	
    v1 = Jrot( v3, -plate[0] )                  # v1 is endpoint for thickness tplate[1]
  path.append( stokes( v1 ) )				    # endpoint of entire stack
  return numpy.array(path)

def makeendpoints( fGHz, vin, stack ) :
  path = [ stokes( vin) ]
  v1 = vin
  for plate in stack :
    v2 = Jrot( v1, plate[0] )                   # v2 is rotated starting point for this segment
    v3 = Jdelay( v2, dphi(L=plate[1], fGHz=fGHz) )      # for sapphire put in different dphi!
    v1 = Jrot( v3, -plate[0] )                  # v1 is endpoint for thickness tplate[1]
    path.append( stokes( v1 ) )				    # endpoint of entire stack
  return numpy.array(path)

def makeequator( ax, viewEl, viewAz ) :
  npts = 256
  x = numpy.zeros( npts )
  y = numpy.zeros( npts )
  z = numpy.zeros( npts )
  n = 0
  for phi in numpy.arange( 0., 2*math.pi, 2*math.pi/npts ) :
    x[n] = math.cos(phi)
    y[n] = math.sin(phi)
    n = n+1
  q = Axes3D.plot(ax, x, y, zs=z, c='gray' )
  #q = Axes3D.plot(ax, x, z, zs=y, c='gray' )
  #q = Axes3D.plot(ax, z, y, zs=x, c='gray' )

  # -- also plot circle that normal to viewing axis -- #
  alpha = math.pi * (90.-viewEl) / 180. 
  xp = x
  yp = math.cos(alpha) * y + math.sin(alpha) * z
  zp = -1.*math.sin(alpha) * y + math.cos(alpha) * z  
  beta = math.pi * (90.-viewAz)/180.
  x = math.cos(beta)*xp + math.sin(beta)*yp 
  y = -1.*math.sin(beta)*xp + math.cos(beta)*yp
  z = zp
  q = Axes3D.plot(ax, x, y, zs=z, c='gray' )

def makeaxes( ax ) :
  x = numpy.array( [-1, 1] )
  y = numpy.array( [0,0] )
  q = Axes3D.plot(ax, x, y, zs=y, c='black' )
  q = Axes3D.plot(ax, y, x, zs=y, c='black' )
  q = Axes3D.plot(ax, y, y, zs=x, c='black' )
  
def addcurve( ax, path, color ) :
  px = path[ :,0 ]
  py = path[ :,1 ]
  pz = path[ :,2 ]
  q = Axes3D.plot(ax, px, py, '--', zs=pz, c=color, linewidth=2 )

def addpoints( ax, endpoints, color ) :
  px = endpoints[ :,0 ]
  py = endpoints[ :,1 ]
  pz = endpoints[ :,2 ]
  q = Axes3D.scatter3D(ax, px,py, zs=pz, c=color, marker='o', s=60, depthshade=False )
     
def doit7( thetadeg ) :
  theta = thetadeg * math.pi/180.
  vin = numpy.array( [math.cos(theta),math.sin(theta)], dtype=complex )
  stack1 = [ [45.,tt] ]
  stack2 = [ [22.5,tt], [67.5,tt] ]
  stack3 = [ [15.,tt],  [75.,tt],   [15.,tt] ]
  stack5 = [ [7.,tt],   [36.,tt],  [102.,tt], [36.,tt], [7.,tt] ]
  which = stack3
  p1 = makepath( 130., vin, which )
  e1 = makeendpoints( 130., vin, which )
  p2 = makepath( 150., vin, which )
  e2 = makeendpoints( 150., vin, which )
  p3 = makepath( 170., vin, which )
  e3 = makeendpoints( 170., vin, which )
  
# shows Faraday rotation
def drawcylinder( ax ) :
  npts = 2048 
  xa = numpy.array( [-1, 1] )
  ya = numpy.array( [0.,0.] )
  za = numpy.array( [40.,40.] )
  q = Axes3D.plot(ax, xa, za, zs=ya, c='gray' )
  q = Axes3D.plot(ax, ya, za, zs=xa, c='gray' )
  za = numpy.array( [-40.,-40.] )
  q = Axes3D.plot(ax, xa, za, zs=ya, c='gray' )
  q = Axes3D.plot(ax, ya, za, zs=xa, c='gray' )
  xa = numpy.array( [0.,0.] )
  za = numpy.array( [0.,0.] )
  ya = numpy.array( [-40.,40.] )
  q = Axes3D.plot(ax, xa, ya, zs=za, c='gray', linestyle='--' )

  x = numpy.zeros( npts )
  y = -30.*numpy.ones( npts )
  z = numpy.zeros( npts )
  xx = numpy.zeros( npts )
  yy = numpy.zeros( npts )
  zz = numpy.zeros( npts )

  n = 0
  for phi in numpy.arange( 0., 2*math.pi, 2*math.pi/npts ) :
    x[n] = math.cos(phi)
    z[n] = math.sin(phi)
    n = n+1
  #q = Axes3D.plot(ax, x, y, zs=z, c='gray' )
  y = 30.* numpy.ones( npts )
  #q = Axes3D.plot(ax, x, y, zs=z, c='gray' )
  n = 0
  for phi in numpy.arange( -30., 30.0, (60./npts) ) :
    yy[n] = phi
    amp = math.cos(phi)
    zz[n] = amp * math.cos( .007 * (30.-phi) )
    xx[n] = -1.* amp * math.sin( .007 * (30.-phi) )
    #zz[n] = math.cos(phi)
    n = n+1
  q = Axes3D.plot(ax, xx, yy, zs=zz, c='black' )
  xxx = numpy.zeros( 2 )
  yyy = 40.*numpy.ones( 2 )
  zzz = numpy.zeros( 2 )
  zzz[0] = -1.
  zzz[1] = 1. 
  q = Axes3D.plot(ax, xxx, yyy, zs=zzz, c='black' )
  xxxx = numpy.zeros( 2 )
  yyyy = -40.*numpy.ones( 2 )
  zzzz = numpy.zeros( 2 )
  xxxx[0] = -1.* math.sin( .007 * (60.) )
  zzzz[0] = math.cos( .007 * (60.) )
  xxxx[1] = 1.* math.sin( .007 * (60.) )
  zzzz[1] = -1.* math.cos( .007 * (60.) )
  q = Axes3D.plot(ax, xxxx, yyyy, zs=zzzz, c='black' )
  

# this shows Faraday rotation!!
def cylinder() :
  fig = plt.figure()
  ax = Axes3D(fig)
  #ax.set_aspect('equal')
  #ax.set_axis_off()
  drawcylinder( ax )
  ax.auto_scale_xyz( [-5,5], [-40,40], [-5,5] )
  ax.set_aspect('equal')
  ax.view_init(25,-75)
  #ax.set_axis_off()
  #q = fig.gca(projection='3d')
  #q._axis3don=False
  plt.show()

# make Poincare plots for 1mm polarizer
# start at theta=0
# specify projections
def doit9( viewEl=8, viewAz=200 ) :
  vin = R
  freqlist = [210.,230.,260.]
  colorlist = ['red','green','blue']
  plt.ion()
  #pp = PdfPages( "Poincare.pdf" )
  stack = [ [-45., length() ] ] 
  #stack = [ [ -15.5, length() ], [ -75., 2.*length() ]  ] 
  fig = plt.figure( figsize=(20,10), facecolor='white', tight_layout=True )
  ax = fig.add_subplot( 1, 2, 1, projection='3d', title="1-section" )
  ax.set_aspect('equal')
  ax.set_axis_off()
  makeequator( ax, viewEl, viewAz )
     # draws equator plus outline of sphere as seen from viewEl and viewAz
     # figure can be rotated on screen, but outline of sphere will then be incorrect,
     # so determine optimum viewEl and viewAz before making final plot
  makeaxes( ax )
     # adds x, y, z axes
  for freq, color in zip( freqlist, colorlist ) :
    p1 = makepath( freq, vin, stack )
    addcurve( ax, p1, color )
  for freq, color in zip( freqlist, colorlist ) :
    e1 = makeendpoints( freq, vin, stack )
    addpoints( ax, e1, color )
  ax.view_init(elev=viewEl,azim=viewAz)

  stack = [ [ -15.5, length() ], [ -75., 2.*length() ]  ] 
  ax = fig.add_subplot( 1, 2, 2, projection='3d', title="2-section" )
  ax.set_aspect('equal')
  ax.set_axis_off()
  makeequator( ax, viewEl, viewAz )
     # draws equator plus outline of sphere as seen from viewEl and viewAz
     # figure can be rotated on screen, but outline of sphere will then be incorrect,
     # so determine optimum viewEl and viewAz before making final plot
  makeaxes( ax )
     # adds x, y, z axes
  for freq, color in zip( freqlist, colorlist ) :
    p1 = makepath( freq, vin, stack )
    #p1 = makepath( 210., vin, [[-75., 1.5*length() ]] )  example of how to do shorter paths
    addcurve( ax, p1, color )
  for freq, color in zip( freqlist, colorlist ) :
    e1 = makeendpoints( freq, vin, stack )
    addpoints( ax, e1, color )
  ax.view_init(elev=viewEl,azim=viewAz)
  #plt.tight_layout()
  plt.show()
  #plt.savefig( pp, format='pdf', bbox_inches='tight' )

# ======================== 4/23/2015 =========================== #

# computes transmission through a birefringent plate for beam incident at angle angIdeg
# this is best done with transfer matrices, but this is the usual FP calculation where the
#        signals keep rattling back and forth
# assume plane of incidence coincides with a principal axis, so that there is no mixing
#        between par and perp polarizations

# ... npar = index of refraction parallel to plane of incidence
# ... nperp = index of refractoin perpendicular to plane of incidence
# ... angIdeg is angle of incidence to the plate
# ... tcm is thickness of plate in cm
# ... alpha is the power absorption coefficient per cm (for both axes)
# ... tanDelta = loss tangent; if non-zero it is used in place of alpha

def plateTrans( fGHz, npar, nperp, angIdeg=45, tcm=1., alpha=0., tanDelta=0. ) :
  print "... npar, nperp = ", npar, nperp
  if tanDelta != 0. :
    alpha = 2. * math.pi * (npar+nperp)/2. * tanDelta * fGHz/clight	  # POWER loss per cm; use avg index for simplicity
  # print "... alpha = %.3f cm^-1" % alpha
  thetaI = math.pi * angIdeg / 180.	                              # convert to radians 
  
# polarization parallel to plane of incidence (n = npar)    
  thetaT = math.asin((math.sin(thetaI))/npar)	                               # Snell's law, eqn 8 (p. 38)
  # print "... angle of incidence = %.2f deg, angle of transmission = %.2f deg" % (angIdeg, 180.*thetaT/math.pi)
  attn = alpha/2. * tcm / math.cos(thetaT) 				
      #  AMPLITUDE attenuation per pass through plate
  [tpar,dummy1,rpar,dummy2] = fresnel( 1., thetaI, npar, thetaT )              # entering dielectric
  [tprimepar,dummy1,rprimepar,dummy2] = fresnel( npar, thetaT, 1., thetaI )    # leaving dielectric
  delta0 = 2 * math.pi * npar * fGHz/clight * tcm / math.cos(thetaT)           # phase shift for 1 pass through plate
  ampTpar = tpar * tprimepar * numpy.exp( -attn + 1.j * delta0 )               # 1st transmitted pass 
  ampRpar = rpar                                                               # 1st reflected ray
  delta2 = 4 * math.pi * npar * fGHz/clight * tcm * math.cos(thetaT)      
      # extra phase shift for each double pass through plate minus phase saved in free space
  for m in range(1,11) :
    # print ".. par  T,R:  %d  %.4f  %.4f" % (m, ampTpar*ampTpar.conjugate(), ampRpar*ampRpar.conjugate() )
    ampTpar = ampTpar + tpar * tprimepar * pow(rprimepar*rprimepar,m) * numpy.exp(-2.*m*attn + 1.j*(delta0 + m*delta2))
    ampRpar = ampRpar + tpar * tprimepar * rprimepar * pow(rprimepar*rprimepar,m-1) * numpy.exp(-2.*m*attn+ 1.j*m*delta2 )
        # each additional pass through the plate

# polarization perpendicular to plane of incidence (n = nperp)
  thetaT = math.asin((math.sin(thetaI))/nperp)	                                         # Snell's law, eqn 8 (p. 38)
  attn = alpha/2. * tcm / math.cos(thetaT)
	# AMPLITUDE attenuation per pass through plate
  [dummy1,tperp,dummy2,rperp] = fresnel( 1., thetaI, nperp, thetaT )              # entering dielectric
  [dummy1,tprimeperp,dummy2,rprimeperp] = fresnel( nperp, thetaT, 1., thetaI )    # leaving dielectric
  # print "... tperp, rperp            ", tperp, rperp
  # print "... tprimeperp, rprimeperp  ", tprimeperp, rprimeperp
  delta0 = 2. * math.pi * nperp * fGHz/clight * tcm / math.cos(thetaT)            # phase shift for 1 pass through plate
  ampTperp = tperp * tprimeperp * numpy.exp( -attn + 1.j * delta0 )               # 1st transmitted pass
  ampRperp = rperp                                                                # 1st reflected pass
  delta2 = 4 * math.pi * nperp * fGHz/clight * tcm * math.cos(thetaT)       
      # extra phase shift from double pass through plate minus phase saved in free space
  print "... delta2 =", 180./math.pi * delta2
  print "... perp phase shift per pass: ", numpy.angle( tprimeperp*rprimeperp * numpy.exp(1.j*delta2), deg=True)
  for m in range(1,11) :
    #print ".. perp T,R:  %d  %.4f  %.4f" % (m, ampTperp*ampTperp.conjugate(), ampRperp*ampRperp.conjugate() )
    ampTperp = ampTperp + tperp * tprimeperp * pow(rprimeperp*rprimeperp,m) * numpy.exp(-2.*m*attn + 1.j*(delta0 + m*delta2))
    ampRperp = ampRperp + tperp * tprimeperp * rprimeperp * pow(rprimeperp*rprimeperp,m-1) * numpy.exp(-2.*m*attn + 1.j*m*delta2 )
  return [ ampTpar, ampTperp, ampRpar, ampRperp ]

def doit2( angIdeg=45., npar=0., nperp=0., fGHz0=220. ) :
  fout = open("plateTrans", "w" )
  if npar==0. :
    npar,nperp = nsapphire( fGHz0 )
  # print "using npar = %.3f, nperp = %.3f" % (npar,nperp)


  fout.write("# angIdeg = %.2f\n" % angIdeg )
  fout.write("# using npar = %.3f, nperp = %.3f\n" % (npar,nperp))
  fout.write("# fGHz     Tpar    Rpar   Apar    Tperp   Rperp   Aperp\n")
  for fGHz in numpy.arange(fGHz0-10.,fGHz0+10.005,.005) :
     ampTpar,ampTperp,ampRpar,ampRperp = plateTrans( fGHz, npar, nperp, angIdeg=angIdeg, alphaO=0.1, alphaE=0.1  )
     Tpar = ampTpar*ampTpar.conjugate()
     Tperp = ampTperp*ampTperp.conjugate()
     Rpar = ampRpar*ampRpar.conjugate()
     Rperp = ampRperp*ampRperp.conjugate()
     fout.write("%.3f    %.4f  %.4f  %.4f    %.4f  %.4f  %.4f\n" % (fGHz, Tpar, Rpar, 1.-Rpar-Tpar, Tperp, Rperp, 1.-Rperp-Tperp ))
     #fout.write("%.3f   %.4f  %.4f  %.4f  %7.2f  %6.4f  %7.2f\n" % (fGHz, ampTpar*ampTpar.conjugate(), ampTperp*ampTperp.conjugate(), \
     #   abs(ampTpar), numpy.angle(ampTpar, deg=True), abs(ampTperp), numpy.angle(ampTperp, deg=True) ) )
  fout.close()

    
def doit3( LOGHz=110., n1=0., n2=0., angIdeg=45., thetaDeg=0., tanDelta=1.e-4, dataFile='junk.dat' ) :
  Thot = 290
  Tcold = 77
  IF = []
  T3 = []
  T4 = []
  T5 = []
  if n1 == 0. : 
    n1 = neff( thetaDeg )
    n2 = neff( 90.+thetaDeg )
  npar = n1
  nperp = n2
  fout = open("ifTrans", "w" )
  fout.write("# ifTrans for LOGHz = %.3f\n" % LOGHz)
  fout.write("# created with pol.doit3()\n")
  fout.write("# plate angle of incidence = %.0f degrees\n" % angIdeg)
  fout.write("# n_parallel = %.2f, n_perpendicular = %.2f\n" % (npar,nperp) )
  fout.write("# IF    trans   refl   loss      T3      T4       T5\n")
  for fIF in numpy.arange(0.2, 10.0, .05 ):
     TparLSB,TperpLSB,RparLSB,RperpLSB = plateTrans( LOGHz-fIF, npar, nperp, angIdeg=angIdeg, tanDelta=tanDelta )
     TparUSB,TperpUSB,RparUSB,RperpUSB = plateTrans( LOGHz+fIF, npar, nperp, angIdeg=angIdeg, tanDelta=tanDelta )
     trans = (TperpLSB*TperpLSB.conjugate() \
               + TperpUSB*TperpUSB.conjugate())/2.
     refl = ( RperpLSB*RperpLSB.conjugate() \
               + RperpUSB*RperpUSB.conjugate())/2.
     #trans = (TparLSB*TparLSB.conjugate() + TperpLSB*TperpLSB.conjugate() \
     #          + TparUSB*TparUSB.conjugate() + TperpUSB*TperpUSB.conjugate())/4.
     #refl = (RparLSB*RparLSB.conjugate() + RperpLSB*RperpLSB.conjugate() \
     #          + RparUSB*RparUSB.conjugate() + RperpUSB*RperpUSB.conjugate())/4.
     loss = 1.-trans-refl
     #print fIF, trans, refl, 1.-trans-refl
     P3 = trans*Thot + refl*Tcold + loss*Thot
     P4 = trans*Tcold + refl*Thot + loss*Thot
     P5 = trans*Tcold + refl*Tcold + loss*Thot
     fout.write("%.3f  %6.4f %6.4f %6.4f   %6.1f  %6.1f  %6.1f\n" % (fIF, trans, refl, loss, P3, P4, P5 ) )
     IF.append( fIF )
     T3.append( P3 )
     T4.append( P4 )
     T5.append( P5 )
  fout.close()
  fin = open( dataFile, "r" )
  fdat = []
  T3dat = []
  T4dat = []
  T5dat = []
  for line in fin :
    if not line.startswith('#') :
      a = line.split()
      fdat.append( float(a[0]) )
      T3dat.append( float(a[1]) )
      T4dat.append( float(a[2]) )
      T5dat.append( float(a[3]) )
  fin.close() 
  plt.ioff()
  fig = plt.subplot(1,1,1)
  fig.plot( IF, T3, color='red' )
  fig.plot( fdat, T3dat, color='red' )
  fig.plot( IF, T4, color='blue' )
  fig.plot( fdat, T4dat, color='blue' )
  fig.plot( IF, T5, color='purple' )
  fig.plot( fdat, T5dat, color='purple' )
  fig.grid( True )
  plt.show()

# compute average power P3,P4,P5 and use it to compute R,T,Abs for a range of tanDelta;
# the idea is to compare the plate absorption with tanDelta

def doit4( nrefrac, LOGHz=220., angIdeg=45., tcm=1. ) :
  Thot = 290
  Tcold = 77
  fout = open("plateLoss", "w" )
  fout.write("# created with pol.doit4()\n")
  fout.write("# LOGHz = %.3f\n" % LOGHz)
  fout.write("# angIdeg = %.2f deg\n" % angIdeg )
  #fout.write("# thetaE = %.2f deg\n" % thetaE )
  fout.write("# nrefrac = %.2f\n" % nrefrac )
  fout.write("# tanDelta     P3     P4     P5       T       R       Abs   tanD_est\n")
  for tanDelta in numpy.arange( 0., 2.05e-3, 2.e-2 ) :
    print "processing tanDelta = %.5f" % tanDelta
    for pp in [ [0.,1.], [1.,0.], [.5,.5] ] :		# [par,perp]
      P3avg = 0.
      P4avg = 0.
      P5avg = 0.
      numpts = 0.
      for fIF in numpy.arange(0.5, 20., .01 ):
        TparLSB,TperpLSB,RparLSB,RperpLSB = plateTrans( LOGHz-fIF, nrefrac, nrefrac, angIdeg=angIdeg, tanDelta=tanDelta, tcm=tcm )
        TparUSB,TperpUSB,RparUSB,RperpUSB = plateTrans( LOGHz+fIF, nrefrac, nrefrac, angIdeg=angIdeg, tanDelta=tanDelta, tcm=tcm )
        trans = (pp[0]*TparLSB*TparLSB.conjugate() + pp[1]*TperpLSB*TperpLSB.conjugate() \
               + pp[0]*TparUSB*TparUSB.conjugate() + pp[1]*TperpUSB*TperpUSB.conjugate())/2.
        refl = (pp[0]*RparLSB*RparLSB.conjugate() + pp[1]*RperpLSB*RperpLSB.conjugate() \
               + pp[0]*RparUSB*RparUSB.conjugate() + pp[1]*RperpUSB*RperpUSB.conjugate())/2.
        loss = 1.-trans-refl
        #print "tanDelta, trans, refl, loss = ", tanDelta, trans, refl, loss    
        P3avg = P3avg + trans*Thot + refl*Tcold + loss*Thot
        P4avg = P4avg + trans*Tcold + refl*Thot + loss*Thot
        P5avg = P5avg + trans*Tcold + refl*Tcold + loss*Thot
        numpts =numpts + 1   
      P3avg = P3avg/numpts
      P4avg = P4avg/numpts
      P5avg = P5avg/numpts
      R = (Thot - P3avg)/(Thot - Tcold)
      T = (Thot - P4avg)/(Thot - Tcold)
      L = (P5avg - Tcold)/(Thot - Tcold)
      thetaT = math.asin((math.sin(angIdeg*math.pi/180.))/nrefrac)	                                         # Snell's law, eqn 8 (p. 38)
      computedAlpha = -1.*math.log(1 - L) * math.cos(thetaT)/tcm
      computedtanDelta = clight*computedAlpha/(LOGHz * 2. * math.pi * nrefrac)
      fout.write("  par %.1f  perp %.1f   %8.5f  %6.1f %6.1f %6.1f  %7.3f %7.3f %7.3f  %7.5f   %7.5f\n" % \
         (pp[0],pp[1],tanDelta, P3avg, P4avg, P5avg, T, R, L, computedAlpha, computedtanDelta) )
  fout.close()

# reads file produced by pb.calcT
def readDataFile( dataFile ) :
  fin = open( dataFile, "r" )
  fdat = []
  T3dat = []
  T4dat = []
  T5dat = []
  for line in fin :
    if not line.startswith('#') :
      a = line.split()
      fdat.append( float(a[0]) )
      T3dat.append( float(a[1]) )
      T4dat.append( float(a[2]) )
      T5dat.append( float(a[3]) ) 
  fin.close() 
  return [ numpy.array(fdat), numpy.array(T3dat), numpy.array(T4dat), numpy.array(T5dat) ]
   
    
def smArray( x, nsm=6 ) :
  smx = []
  i = 0
  while (i + nsm) < len(x) :
    smx.append( numpy.mean( x[i:i+nsm] ) )
    i = i + nsm
  return( numpy.array( smx ) )
    
# doit5 tries to FIT the data
# angIrange is angle of incidence of the plate (nominally 45 degrees)
# thetaRange is angle of ordinary axis relative to plane of incidence, in the frame of the plate

def doit5( dataFile, angIrange=[43.,48.,2.], phiRange=[-5.,6.,2.], alphaO=.05, alphaE=.03 ) :
  Thot = 295
  Tcold = 77
  varsave = 1.e6
  LOGHz = float( dataFile.split("_")[0] )
  print "LOGHz = %0f" % LOGHz
  fdat, T3dat, T4dat, T5dat = readDataFile( dataFile )

  # to save time in fitting, boxcar average to shorter arrays
  fsm = smArray( fdat )
  # print fsm    # print freqs in smooth array
  T3sm = smArray( T3dat )
  T4sm = smArray( T4dat )
  T5sm = smArray( T5dat )
  for phi in numpy.arange( phiRange[0], phiRange[1], phiRange[2] ) :
    for angI in numpy.arange( angIrange[0], angIrange[1], angIrange[2] ) :
      #T3,T4,T5 = Tfit( LOGHz, fsm, angI, phi, alphaO, alphaE )  
      T3,T4,T5 = Tfit2( LOGHz, fsm, angI, phi, alphaO, alphaE )  
      variance = numpy.mean( pow( T3-T3sm, 2. )) + \
         numpy.mean( pow( T4-T4sm, 2. ))
      print "phi = %.1f, angI = %.1f, variance = %.1f" % (phi, angI, variance)
      if variance < varsave :
        angIbest = angI
        phibest = phi 
        varsave = variance

  # once best fit is found, model all the points
  var5 = numpy.mean( pow(T5-T5sm, 2.) )
  print "lowest variance = %.1f, var5 = %.2f, for phi  = %.1f, angI = %.1f" % (varsave, var5, phibest, angIbest)
  T3,T4,T5 = Tfit2( LOGHz, fdat, angIbest, phibest, alphaO, alphaE )  
  plt.ioff()
  fig = plt.subplot(1,1,1)
  fig.plot( fdat, T3, color='red' )
  fig.plot( fdat, T3dat, color='red' )
  #fig.plot( fsm, T3sm, color='red' )
  fig.plot( fdat, T4, color='blue' )
  fig.plot( fdat, T4dat, color='blue' )
  #fig.plot( fsm, T4sm, color='blue' )
  fig.plot( fdat, T5, color='purple' )
  fig.plot( fdat, T5dat, color='purple' )
  #fig.plot( fsm, T5sm, color='purple' )
  fig.grid( True )
  plt.title( dataFile )
  fig.text( .05, .05, "phi = %.1f, angI = %.1f, alphaO = %.3f, alphaE = %.3f" % (phibest,angIbest,alphaO,alphaE), \
     transform=fig.transAxes, horizontalalignment="left", verticalalignment="center", size="medium" )
  plt.show( )

# Tfit is the original fitting routine
# it computes npar and nperp from phiDeg, but does NOT keep track of polarization state within 
#	 the plate - that is, it assumes that input polarization is it

def Tfit( LOGHz, fdat, angIdeg, phiDeg, alphaO, alphaE ) :
  Thot = 295.
  Tcold = 77.
  npar = neff( phiDeg )
  nperp = neff( phiDeg + 90. )
  tanDelta=4.e-4
  print "using npar = %.3f, nperp = %.3f" % (npar,nperp)
  T3 = []
  T4 = []
  T5 = []
  
  for fIF in fdat : 
    TparLSB,TperpLSB,RparLSB,RperpLSB = plateTrans( LOGHz-fIF, npar, nperp, angIdeg=angIdeg, tanDelta=tanDelta )
    TparUSB,TperpUSB,RparUSB,RperpUSB = plateTrans( LOGHz+fIF, npar, nperp, angIdeg=angIdeg, tanDelta=tanDelta )
    if ( LOGHz < 200. ) :    # single polarization perp to plane of incidence
      trans = (TperpLSB*TperpLSB.conjugate() \
         + TperpUSB*TperpUSB.conjugate())/2.
      refl = ( RperpLSB*RperpLSB.conjugate() \
         + RperpUSB*RperpUSB.conjugate())/2.
    else :                   # circular pol
      trans = (TparLSB*TparLSB.conjugate() + TperpLSB*TperpLSB.conjugate() \
         + TparUSB*TparUSB.conjugate() + TperpUSB*TperpUSB.conjugate())/4.
      refl = (RparLSB*RparLSB.conjugate() + RperpLSB*RperpLSB.conjugate() \
         + RparUSB*RparUSB.conjugate() + RperpUSB*RperpUSB.conjugate())/4.
    loss = 1.-trans-refl
    # print "... %.3f GHz: T,R,L = %.4f, %.4f, %.4f" % (fIF, trans, refl, loss)
    T3.append( trans*Thot + refl*Tcold + loss*Thot )
    T4.append( trans*Tcold + refl*Thot + loss*Thot )
    T5.append( trans*Tcold + refl*Tcold + loss*Thot )
  return [ numpy.array(T3,dtype=float), numpy.array(T4,dtype=float), numpy.array(T5,dtype=float ) ]

# Tfit2 is the new fitting routine (uses plateTrans2)
# it rotates into the frame of the principal axes inside the plate, keeps track of polarization

def Tfit2( LOGHz, fdat, angIdeg, phiDeg, alphaO, alphaE ) :
  Thot = 295.
  Tcold = 77.
  T3 = []
  T4 = []
  T5 = []
  inpolList = [ Y ]	        	# 3mm band is sensitive to single pol in Y-direction; necessary!!
  if LOGHz > 180. :
    inpolList = [ X, Y ]        # 1mm band is sensitive to circular pol

  for fIF in fdat :             
    trans = 0.
    refl = 0.
    for fGHz in [ LOGHz-fIF, LOGHz+fIF ] :	        # process each sideband separately
      for inpol in inpolList : 					    # process each polarization separately
        Tpar,Tperp,Rpar,Rperp = plateTrans2( inpol, fGHz, angIdeg, phiDeg, alphaO, alphaE )
        trans = trans + Tpar*Tpar.conjugate() + Tperp*Tperp.conjugate()
          # compute power transmitted through plate in both X and Y polarizations
        refl = refl + Rpar*Rpar.conjugate() + Rperp*Rperp.conjugate()
          # compute power reflected from plate in both X and Y polarizations
    for j in range(0,len(inpolList)) :      # divide by 2 or by 4
      trans = trans/2. 
      refl = refl/2.
    loss = 1.-trans-refl
    # print "... %.3f GHz: T,R,L = %.4f, %.4f, %.4f" % (fIF, trans, refl, loss)
    T3.append( trans*Thot + refl*Tcold + loss*Thot )
    T4.append( trans*Tcold + refl*Thot + loss*Thot )
    T5.append( trans*Tcold + refl*Tcold + loss*Thot )
  return [ numpy.array(T3,dtype=float), numpy.array(T4,dtype=float), numpy.array(T5,dtype=float ) ]


# compute effective indices of refraction parallel and perpendicular to plane of incidence
# for a uniaxial dielectric slab in X'-Y-Z' coordinate frame aligned with the slab
# ... phiDeg is angle, in the frame of the plate, of C-axis relative to plane of incidence;
#       if phiDeg=0, nperp=nO; npar ~ nE ; par is NOT parallel with C-axis, hence n is NOT
#       exactly equal to nE, though it is close; if phiDeg=90, nperp=nE, npar = nO
# ... thetaDeg is angle of propagation of the ray in the X-Z plane, relative to Z-axis,
#       INSIDE the dielectric

def effindex( phiDeg, thetaDeg, nO=3.05, nE=3.40 ) :
  theta = math.pi/180. * thetaDeg
  phi = math.pi/180. * phiDeg
  #vecOrdinary = numpy.array( [ math.cos(phi), math.sin(phi), 0. ], dtype=float )
	# this is unit vector in direction of the ordinary axis
  #vecRay = numpy.array( [ 0., math.sin(theta), math.cos(theta) ] )
    # this is unit vector in the direction of the ray
  #cosA = numpy.dot( vecOrdinary, vecRay )  # law of cosines
    # this is the angle between the ray direction and the ordinary axis
  cosphi = math.cos( phi )
  sinphi = math.sin( phi )
  npar = 1./ math.sqrt( pow(sinphi/nO, 2.) + pow( cosphi/nE, 2. ) )
  nperp = 1./ math.sqrt( pow(cosphi/nO, 2.) + pow( sinphi/nE, 2.) ) 
  return npar,nperp
  
# propagate signal through anisotropic dielectric, at angle phiDeg to ordinary axis
# break wave into components along principal axes, demonstrate that there is an effective
#   index of refraction for this wave
def dumbtest( phiDeg, no=1.5, ne=3.) :
  fout = open( "dumbtest", "w" )
  v1 = X  
  v2 = Jrot( v1, phiDeg )
  for delay1Deg in numpy.arange(0.,3602.,2.) :
    delay2Deg = delay1Deg/1.1
    v3 = J3delay( v2, delay1Deg, delay2Deg )
    v4 = Jrot( v3, -1.*phiDeg)
    fout.write("%.2f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f %8.3f\n" % (delay1Deg, abs(v3[0]), numpy.angle(v3[0], deg=True), \
       abs(v3[1]), numpy.angle(v3[1],deg=True), abs(v4[0]), numpy.angle(v4[0],deg=True)))
  fout.close()

def dumbtest2( no=3.05, ne=3.40 ) :
  fout = open( "dumbtest2", "w" )
  v1 = X  
  tcm = 28.
  for phiDeg in numpy.arange(0.,91.,1.) :
    print "\nv1: ",v1
    v2 = Jrot( v1, phiDeg )
    print "v2: ",v2
    delay1Deg = 360.*220.e9 * no * tcm/3.e10   
    #delay1Deg = 360090.
    print "delay1Deg =", delay1Deg
    delay2Deg = 360.*220.e9 * ne * tcm/3.e10
    #delay2Deg = 360180.
    print "delay2Deg =", delay2Deg
    v3 = J3delay( v2, delay1Deg, delay2Deg )
    print "v3: ",v3
    v4 = Jrot( v3, -1.*phiDeg)
    print "v4: ",v4
    phi = phiDeg * math.pi/180.
    sq = pow( math.cos(phi)/no,2.) + pow(math.sin(phi)/ne,2.) 
    neff = 1./math.sqrt(sq)
    #phpredict = numpy.mod( 360.*220.e9 * neff * 1./3.e10, 360.)
    phpredict =  360.*220.e9 * neff * tcm/3.e10
    fout.write("%.1f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n" % \
       (phiDeg, abs(v4[0]), numpy.angle(v4[0], deg=True), neff, phpredict, delay1Deg, delay2Deg))
  fout.close()

def dumbtest3( ) :
  fout = open("dumbtest3", "w")
  x = numpy.arange(0.,50.,.01)
  v1 = numpy.cos( 3*x )
  v2 = numpy.cos( 2*x )
  v3 = v1 + v2
  v4 = numpy.cos( 2.5 * x)
  v5 = numpy.cos( 2.353394 * x)
  for xx,vv1,vv2,vv3,vv4,vv5 in zip(x,v1,v2,v3,v4,v5) :
    fout.write( "%8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n" % (xx,vv1,vv2,vv3,vv4,vv5) )
  fout.close

# compute amplitude transmission and reflection through surface
# the input Jones vector v1 = [ Ex, Ey] = [ Epar, Eperp ] where par and perp are relative to the plane of incidence
# return trans and refl Jones vectors
def surf( vec, tpar, tperp, rpar, rperp ) :
  Tmat = numpy.array( [[ tpar, 0.],[ 0., tperp]] )
  Rmat = numpy.array( [[ rpar, 0.],[ 0., rperp]] )
  return [ numpy.dot(Tmat,vec), numpy.dot(Rmat,vec) ]
  

# compute transmission through uniaxial birefringent plate
# neither principal axis need be aligned with plane of incidence
# ... inpol is the Jones vector [Ex Ey] of incident electric field, in X-Y-Z frame aligned 
#       with X-Z plane of incidence; thus Ey is polarized perp to plane of incidence
# ... angIdeg is the angle of incidence onto the plate
# ... tcm is thickness of plate in cm
# ... phiDeg is angle, in the frame of the plate, of C-axis relative to plane of incidence;
#       if phiDeg=0, nperp=nO; npar ~ nE ; par is NOT parallel with C-axis, hence n is NOT
#       exactly equal to nE, though it is close; if phiDeg=90, nperp=nE, npar = nO
# ... nO = ordinary index of refrac (default = 3.05 for sapphire)
# ... nE = extraordinary index of refrac (default = 3.40 for sapphire)

# here's the plan:
#  - transmission, reflection from surfaces is done in X-Y-Z coordinate system aligned with
#      the plane of incidence, using effective indices of refrac computed from index ellipsoid
#  - calculation of phase shifts and loss through the plate are done in the frame of the
#      principal axes of the plate

def plateTrans2( inpol, fGHz, angIdeg=45., phiDeg=0., alphaO=0., alphaE=0.,  nO=3.05, nE=3.40, tcm=1. ) :
  angI = math.pi * angIdeg / 180.	                              # convert to radians 
  neffpar,neffperp = effindex( phiDeg, angIdeg, nO=nO, nE=nE )
  print "... phiDeg = %.1f  npar,nperp = %.3f, %.3f" % (phiDeg,neffpar,neffperp)

  # compute Fresnel coefficients for polarization parallel to plane of incidence     
  angTpar = math.asin((math.sin(angI))/neffpar)	                          # Snell's law, eqn 8 (p. 38)
  [tpar,dummy1,rpar,dummy2] = fresnel( 1., angI, neffpar, angTpar )              # entering dielectric
  [tprimepar,dummy1,rprimepar,dummy2] = fresnel( neffpar, angTpar, 1., angI )    # leaving dielectric
  
  # compute Fresnel coefficients for polarization perpendicular to plane of incidence
  angTperp = math.asin((math.sin(angI))/neffperp)                             # Snell's law, eqn 8 (p. 38)
  [dummy1,tperp,dummy2,rperp] = fresnel( 1., angI, neffperp, angTperp )                 # entering dielectric
  [dummy1,tprimeperp,dummy2,rprimeperp] = fresnel( neffperp, angTperp, 1., angI )       # leaving dielectric
  avAngT = math.sqrt( angTpar*angTperp )

  # compute atten and phase shifts (in radians) through the plate
  angTO = math.asin((math.sin(angI))/nO)
  angTE = math.asin((math.sin(angI))/nE)
  #attnO = alphaO/2. * tcm / math.cos(angTO) 					#  AMPLITUDE attenuation per pass through plate
  #attnE = alphaE/2. * tcm / math.cos(angTE) 					#  AMPLITUDE attenuation per pass through plate
  #deltaO = 2 * math.pi * nO * fGHz/clight * tcm / math.cos(angTO)     
  #deltaE = 2 * math.pi * nE * fGHz/clight * tcm / math.cos(angTE)    
  attnO = alphaO/2. * tcm / math.cos(avAngT) 					#  AMPLITUDE attenuation per pass through plate
  attnE = alphaE/2. * tcm / math.cos(avAngT) 					#  AMPLITUDE attenuation per pass through plate
  deltaO = 2 * math.pi * nO * fGHz/clight * tcm / math.cos(avAngT)     
  deltaE = 2 * math.pi * nE * fGHz/clight * tcm / math.cos(avAngT)    
  # print "... deltaO, deltaE = ", 180./math.pi * deltaO, 180./math.pi * deltaE
  delaymat = numpy.array( [ [cmath.exp(-attnE - 1j*deltaE), 0], [0, cmath.exp(-attnO - 1j*deltaO)] ] )
  # extraphasePerp = math.sin(angI) * 2. * tcm * math.tan(angTperp) * 2. * math.pi * fGHz/clight
  # extraphasePar = math.sin(angI) * 2. * tcm * math.tan(angTpar) * 2. * math.pi * fGHz/clight
  extraphasePerp = math.sin(angI) * 2. * tcm * math.tan(avAngT) * 2. * math.pi * fGHz/clight
  extraphasePar = math.sin(angI) * 2. * tcm * math.tan(avAngT) * 2. * math.pi * fGHz/clight
  print "... extraphasePar, Perp = ", 180./math.pi * extraphasePar, 180./math.pi * extraphasePerp
    # this is the phase SAVED outside the plate on each pass due to geometry; depends on neffpar and neffperp
  extramat = numpy.array( [ [cmath.exp(1j*extraphasePar), 0], [0, cmath.exp(1j*extraphasePerp)] ] )
  transCorr = numpy.array( [ [neffpar * math.cos(angTpar)/math.cos(angI), 0.], [0., neffperp * math.cos(angTperp)/math.cos(angI)] ] )
    # this is the irradiance correction factor just inside the surface, in the X-Y frame only

  # begin with Jones vector in frame of tilt axis
  v1 = inpol 
  v2,R = surf( v1, tpar, tperp, rpar, rperp )                # v2 is Jones vector just inside plate, R is initial reflected amp
  T = numpy.array( [0.,0.], dtype=complex )                  # Jones vector transmitted; each component is complex
  Acalc = numpy.array( [0.,0.], dtype=float )
  Rsave = R

  # now let signal rattle back and forth inside the plate; each pass is one round trip
  for npass in range(0, 12) :  
    v3 = Jrot( v2, phiDeg )                     # rotate into frame of principal axes [ E_O, E_E ]
    v4 = numpy.dot( delaymat, v3 )              # complex delays include attenuation
    v5 = Jrot( v4, -phiDeg )                    # rotate back to original X-Y frame
    Acalc = Acalc + numpy.dot( transCorr, v2*v2.conjugate() - v5*v5.conjugate() )
    t,v2 = surf( v5, tprimepar, tprimeperp, rprimepar, rprimeperp)    
       # t is the Jones vector transmitted through the surface
       # v2 is the Jones vector reflected just inside the surface
    extramat = numpy.array( [ [cmath.exp(1j*npass*extraphasePar), 0], [0, cmath.exp(1j*npass*extraphasePerp)] ] )
    t = numpy.dot( extramat, t )     
       # apply geometric phase shift to both components of t; no correction for pass0 = reference phase
    T = T + t                                   # add to T (in X,Y frame)
    v3 = Jrot( v2, phiDeg )                     # direction of propagation is reversed; should I reverse phi ??
    v4 = numpy.dot( delaymat, v3 )              # complex delays include attenuation
    v5 = Jrot( v4, -phiDeg )                    # rotate back to original X-Y frame
    Acalc = Acalc + numpy.dot( transCorr, v2*v2.conjugate() - v5*v5.conjugate() )
    t,v2 = surf( v5, tprimepar, tprimeperp, rprimepar, rprimeperp)    
       # compute transmission, reflection from 1st surface; note that transmitted wave will add to reflection
    extramat = numpy.array( [ [cmath.exp(1j*(npass+1)*extraphasePar), 0], [0, cmath.exp(1j*(npass+1)*extraphasePerp)] ] )
    t = numpy.dot( extramat, t)
       # (npass+1)*extraphase is geometric delay outside the plate SAVED after npasses
    R = R + t 								    # at the first surface, transmitted signals add to R

    # --- debug info: make sure that all power is accounted for ---
    T0tmp = T[0]*T[0].conjugate()
    T1tmp = T[1]*T[1].conjugate()
    R0tmp = R[0]*R[0].conjugate()
    R1tmp = R[1]*R[1].conjugate()
    Atmp = 1. - T0tmp - T1tmp - R0tmp - R1tmp 
    print "... after %2d passes, T = %.5f, %.5f, R = %.5f, %.5f, A = %8.5f" % \
       (npass, abs(T0tmp),abs(T1tmp),abs(R0tmp),abs(R1tmp), Atmp )
    # print "... after %2d passes, T = %.5f, %.5f, R = %.5f, %.5f, Acalc = %.5f, %.5f, A = %8.5f" % \
     #   (npass, abs(T0tmp),abs(T1tmp),abs(R0tmp),abs(R1tmp), Acalc[0], Acalc[1], Atmp )
  return [ T[0], T[1], R[0], R[1] ]             # [ ampTpar, ampTperp, ampRpar, ampRperp ]

# compares results of plateTrans and plateTrans2 for special case where phiDeg=0 or phiDeg=90 only

def ptcompare( LOGHz, angIdeg, phiDeg0=True, alpha=0. ) :
  v1 = numpy.array( [0,1], dtype=complex )
  nO = 3.05
  nE = 3.40
  npar = nE
  nperp = nO
  phiDeg = 0.
  if not phiDeg0 :
    phiDeg = 90. 
    npar = nO
    nperp = nE
  print "\nplateTrans:"
  Tpar,Tperp,Rpar,Rperp = plateTrans( LOGHz, npar, nperp, angIdeg=angIdeg, alpha=alpha )
  #print "%8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f   %5.3f %5.3f %5.3f %5.3f" % (Tpar.real,Tpar.imag, \
  #  Tperp.real,Tperp.imag,Rpar.real,Rpar.imag,Rperp.real,Rperp.imag, abs(Tpar), abs(Tperp), abs(Rpar), abs(Rperp) )
  T = abs(Tperp*Tperp.conjugate()) 
  R = abs(Rperp*Rperp.conjugate())
  #T = abs(Tpar*Tpar.conjugate() + Tperp*Tperp.conjugate()) / 2.
  #R = abs(Rpar*Rpar.conjugate() + Rperp*Rperp.conjugate()) / 2.
  A = 1 - T - R
  print "... T = %.4f, R = %.4f, A = %.4f" % (T, R, A)
  print "\nplateTrans2:"
  Tpar,Tperp,Rpar,Rperp = plateTrans2( v1, LOGHz, angIdeg, phiDeg, alpha, alpha, nO, nE )
  #print "%8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f   %8.4f %8.4f   %5.3f %5.3f %5.3f %5.3f" % (Tpar.real,Tpar.imag, \
  #  Tperp.real,Tperp.imag,Rpar.real,Rpar.imag,Rperp.real,Rperp.imag, abs(Tpar), abs(Tperp), abs(Rpar), abs(Rperp) )
  T = abs(Tpar*Tpar.conjugate() + Tperp*Tperp.conjugate())
  R = abs(Rpar*Rpar.conjugate() + Rperp*Rperp.conjugate())
  A = 1 - T - R
  print "... T = %.4f, R = %.4f, A = %.4f" % (T, R, A)
