# ---------------------------------------------------------------------------------------------------------- #
# pol.py

#from Numeric import *
from numpy import *
import math
import cmath
import random

r_default = 0.061
f_default = 0.015
Rc_default = 0.07
Ltaper_default = 0.08

# a polarizer is a LIST of retarder sections, each of which is a DICTIONARY
# section 1 of the polarizer is closest to the OMT
#
# retarder dictionary keywords: 
#   "angle" = total rotation relative to OMT Y-axis, viewed from OMT toward feed horn
#   "r" = radius of circular waveguide
#   "f" = facet depth of flats
#   "Ltaper" = length of linear taper at each end of the flat, IF Rc WERE 0!
#   "Lflat" = length of the flat section, not including tapers, IF Rc WERE 0!
#   "Rc" = fillet radius at intersection of taper and flat 
#
# the following are DERIVED quantities:
#   "LtaperLin" = actual length of linear taper, not including fillet
#   "Lfillet" = length of filleted section
#   "LflatTrue" = actual length of flat section, not including fillet

# ---------------------------------------------------------------------------------------------------------- #
# define basis vectors in reference coordinate system (aligned with OMT axes)
# each vector consists of 2 complex numbers [ (Re_x,Im_x), (Re_y,Im_y) ]
# R and L propagate in +Z direction according to right hand rule
# note: x.conjugate() gives complex conjugate of a number
# ---------------------------------------------------------------------------------------------------------- #
X = array( [1,0] )
Y = array( [0,1] )
R = array( [1./sqrt(2.), -1j*(1./sqrt(2.))] ) 
L = array( [1./sqrt(2.), +1j*(1./sqrt(2.))] ) 
clight = 29.9792458	# speed of light, cm/nanosec

# ---------------------------------------------------------------------------------------------------------- #
# returns new basis vector in a coordinate system rotated by thetaDegrees
# ---------------------------------------------------------------------------------------------------------- #
def Jrot ( vec, thetaDegrees ) :
  rad = math.pi * thetaDegrees/180.
  rotmat = array( [[cos(rad),sin(rad)],[-sin(rad),cos(rad)]] )
  return dot(rotmat,vec)

# ---------------------------------------------------------------------------------------------------------- #
# returns basis vector after passing through polarizer section
# component 2 (Y-axis) is advanced by delayDegrees relative to component 1 (X-axis)
# ---------------------------------------------------------------------------------------------------------- #
def Jdelay ( vec, delayDegrees ) :
  rad = math.pi * delayDegrees/180.
  rotmat = array( [ [1, 0], [0, cmath.exp(-1.j * rad)] ] )
  return dot( rotmat,vec )

# ---------------------------------------------------------------------------------------------------------- #
# returns basis vector after transmission through a beamsplitter
# tpel is thickness in inches, angI is angle of incidence in degrees
# X axis is assumed parallel to the plane of incidence; light polarized parallel to the plane of incidence
#   is less strongly reflected; at Brewster's angle none will be reflected
# Y axis is perp to plane of incidence; it is more strongly reflected
# ---------------------------------------------------------------------------------------------------------- #
def Jbsplit ( vec, tpel=.001, angI=45., fGHz=90. ) :
   [tpar,tperp,Rpar,Rperp] = pellicle( tpel=tpel, angI=angI, fGHz=fGHz)
   rotmat = array( [ [tpar, 0.], [0., tperp] ] )
   return dot( rotmat,vec ) 

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
def length( r=r_default, f=f_default, dphi0=90., fGHz=85., fcX=0., fcY=0. ) :
  if (f == 0.) : return 0.
  if ( (fcX == 0.) or (fcY == 0.) ) :
    [fcX,fcY] = cutoff( r, f/r )
  length = clight * dphi0 / (2.54 * 360.* (sqrt(fGHz*fGHz-fcY*fcY) - sqrt(fGHz*fGHz-fcX*fcX)))
  return length

# ---------------------------------------------------------------------------------------------------------- #
# compute phase delay in degrees of Y-pol vs X-pol signals traveling through faceted section of length L,
#   radius r, facet depth f; dimensions in inches
# fcX and fcY are override cutoff freqs in GHz; if either is zero, use values from subroutine cutoff
# ---------------------------------------------------------------------------------------------------------- #
def dphi( r=r_default, f=f_default, L=0.001, fGHz=85., fcX=0., fcY=0. ) :
  if ( (fcX == 0.) or (fcY == 0.) ) :
    [fcX,fcY] = cutoff( r, f/r )
  dphi = 360. * L * 2.54 * (sqrt(fGHz*fGHz-fcY*fcY) - sqrt(fGHz*fGHz-fcX*fcX)) / clight
  return dphi

# ---------------------------------------------------------------------------------------------------------- #
# compute differential phase shift through linear tapered section at frequency fGHz
# the facet depth starts at 0.000 and increases linearly to f; all dimensions in inches
# ---------------------------------------------------------------------------------------------------------- #
def dphiTaper( r=r_default, f=f_default, Ltaper=Ltaper_default, fGHz=85., nsteps=100 ) :
  if Ltaper == 0. :
    return 0.
  dL = Ltaper/nsteps      # increment in x
  df = f/nsteps      # increment in facet depth
  x = dL/2.          # mean x-coordinate of 1st section
  fx = f - df/2.     # mean facet depth of 1st section
  totphase = 0.
  for n in range(nsteps) :
    totphase = totphase + dphi(r=r, f=fx, L=dL, fGHz=fGHz)
    fx = fx - df
  return totphase
  #return (totphase + 1.)
    # 15jan - added 1 degree to get closer to HFSS simulations!

# ---------------------------------------------------------------------------------------------------------- #
# compute [length, differential phase shift] through a filleted section
# Rc is the cutter radius in inches; axis of cutter is perpendicular to axis of waveguide
# starting depth is fstart, stopping depth is fstop (use 0 if fillet extends to outer diameter)
# ---------------------------------------------------------------------------------------------------------- #
def dphiFillet ( r=r_default, fstart=f_default, fstop=0., Rc=0.0625, fGHz=85., nsteps=100 ) :
  if Rc == 0. :
    return 0.
  if (fstart-fstop) > Rc :
    fstop = fstart - Rc
    print "WARNING: resetting fstop to fstart - Rc"              # protect against depth > Rc
  totphase = 0.
  df = (fstart - fstop)/nsteps                                   # increment in facet depth
  f1 = fstart                                                    # facet depth at start of segment
  x1 = 0.                                                        # x-coordinate at start of segment
  for n in range(nsteps) :
    f2 = f1 - df                                                 # facet depth at end of segment
    x2 = math.sqrt(Rc*Rc - math.pow((Rc-(fstart-f2)), 2.))       # x-coordinate at end of segment
    dphase = dphi(r=r, f=(f1+f2)/2., L=(x2-x1), fGHz=fGHz )      # using avg facet depth
    totphase = totphase + dphase
    x1 = x2
    f1 = f2
  #return [x2, totphase]           # return X length of transition in inches, and total diff phase through it
  return totphase

# ---------------------------------------------------------------------------------------------------------- #
# reflection coefficient from circular guide of radius r1 (inches) to guide of radius r2 (inches)
# ---------------------------------------------------------------------------------------------------------- #
def reflect( fGHz, r1, r2 ) :
  fc1 = clight/(3.4126 * r1 * 2.54)
  fc2 = clight/(3.4126 * r2 * 2.54)
  Z1 = fGHz * 377./sqrt(fGHz*fGHz - fc1*fc1)
  Z2 = fGHz * 377./sqrt(fGHz*fGHz - fc2*fc2)
  vR = (Z1 - Z2)/(Z1 + Z2)
  print "voltage, power reflection coefficients: %.3e %.3e" % (vR, vR*vR)

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



def blindSearch() :

  r = r_default
  f = f_default
  L90 = length( r=r, f=f, dphi0=90., fGHz=90. ) 
  print "L90 = ", L90
  L180 = 2. * L90
  bestleak = 1.

  #a1 = 6.50
  #a2 = 28.07
  #a3 = 66.57

  for a1 in arange(5.0, 7.05, .05) :
    for a2 in arange( 27.0, 29.05, .05) :
      for a3 in arange (64.0, 68.05, .05) :
        dumpDimensions( 'dim.dat', r, a1, f, L180, a2, f, L180, a3, f, L90 )
        maxleak = computeLeakage( 'dim.dat', 'leak.dat' )
        if (maxleak < bestleak) :
          bestleak = maxleak
          print "best so far: a1 = %.2f, a2 = %2f, a3 = %.2f; bestleak = %.3f" % (a1,a2,a3,bestleak)
   

# ---------------------------------------------------------------------------------------------------------- #

def blindSearch2() :

  r = r_default
  f = f_default
  L90 = length( r=r, f=f, dphi0=90., fGHz=85. ) 
  L180 = 2. * L90
  bestleak = 1.

  L1 = L180
  L2 = L180  
  L3 = L90

  #a1 = 6.50
  #a2 = 28.07
  #a3 = 66.57
 
  #a1 = 6.95
  #a2 = 27.4
  #a3 = 65.45
  
  a10 = 7.5
  a20 = 26.8
  a30 = 64.3

  for a1 in arange( a10-.10, a10+.11, .05) :
    for a2 in arange( a20-.10, a20+.11, .05) :
      for a3 in arange( a30-.10, a30+.11, .05) :
        for L1 in arange( L180-.01, L180+.011, .001) :
          for L2 in arange( L180-.01, L180+.011, .001) :
             for L3 in arange( L90 -.010, L90+.011, .001) :
                maxleak = 0.
                for fGHz in arange( 75., 116., 5.) :
                   leak = quickLeak( r, a1, f, L1, a2, f, L2, a3, f, L3, fGHz )
                   # print fGHz, leak, bestleak
                   if (leak > maxleak) :
                     maxleak = leak
                   if maxleak > bestleak :
                     # print "exiting loop; leak > bestleak"
                     break
                if (maxleak < bestleak) :
                  bestleak = maxleak 
                  # print "setting bestleak = ", bestleak
                  print "best so far: %.2f %3f %2f %3f %.2f %3f  %.4f" % \
                    (a1,L1,a2,L2,a3,L3,bestleak)

def blindSearch3() :

  r = r_default
  f = f_default
  L90 = length( r=r, f=f, dphi0=90., fGHz=85. ) 
  L180 = 2. * L90
  bestleak = 1.

  L1 = L180
  L2 = L180  
  L3 = L90

  #a1 = 6.50
  #a2 = 28.07
  #a3 = 66.57
 
  #a1 = 6.95
  #a2 = 27.4
  #a3 = 65.45
  
  a10 = 7.5
  a20 = 26.8
  a30 = 64.3

  # best so far: 1.00 0.100000 1.000000 0.270000 45.00 0.190000  0.2829

  maxallowedleak = 0.10

  fout = open("blindsearch5.dat", "a")
  fout.write("# restarting search\n")
  fout.flush()
  for a1 in arange( 55., 57., 2. ) :   # i.e., only 55
    #a2targ = 78. - 0.38*a1
    a2targ = 43.
    for a2 in arange( a2targ-1., a2targ+1.,0.2 ) :
      for a3 in arange( 106.8, 107.4, 0.2 ) :
        for L1 in arange( .10, .50, .015 ) :
          for L2 in arange( .10, .50, .015 ) :
             for L3 in arange( .10, .50, .015 ) :
                maxleak = 0.
                for fGHz in arange( 75., 116., 5.) :
                   leak = quickLeak( r, a1, f, L1, a2, f, L2, a3, f, L3, fGHz )
                   # print fGHz, leak, bestleak
                   if (leak > maxleak) :
                     maxleak = leak
                   if (maxleak > bestleak) and (maxleak > maxallowedleak) :
                     # print "exiting loop; leak > bestleak"
                     break
                if (maxleak < bestleak) :
                  bestleak = maxleak 
                  # print "setting bestleak = ", bestleak
                  print "best so far: %.2f %3f %2f %3f %.2f %3f  %.4f" % \
                    (a1,L1,a2,L2,a3,L3,bestleak)
                if (maxleak < maxallowedleak ) :
                  fout.write( "%6.2f %6.3f %6.2f %6.3f %6.2f %6.3f   %6.4f\n" % \
                    (a1,L1,a2,L2,a3,L3,maxleak) )
                  fout.flush()
      fout.write("# finished all lengths, all a3, for a1,a2 = %.2f %.2f\n" % (a1,a2)) 
      fout.flush()
  fout.close()

# -----------------------------------------------------------------------------------

# creates dictionary for a single polarizer section:
#    "r"      is the circular waveguide radius (inches)
#    "f"      is the facet depth (inches)
#    "ang"    is the angle relative to the Yaxis of OMT, viewed from OMT 
#    "Lflat"  is the nominal length of flat (inches), if Rc were 0
#    "Ltaper" is the nominal length of the linear taper (inches), if Rc were 0
#    "Rc"     is the fillet radius between taper and flat
# other quantities that can be derived from these:
#    theta = arctan(f/Ltaper)
#    Lfillet =  length of each fillet section = Rc * math.sin(theta)
#    LtaperLin = length of linear taper = Ltaper * Rc * (1-cos(theta)) / f
#    LflatTrue = length of true flat = Lflat - 2. * (Rc * sin(theta) - Ltaper * Rc (1.-cos(theta)) /f )

def oneSection( ang=0., Lflat=0., LflatTrue=0., r=r_default, f=f_default, Ltaper=Ltaper_default, Rc=Rc_default ) :
  theta = math.atan2(f,Ltaper)
      # slope angle of the taper
  Lfillet = Rc * math.sin(theta)
      # length of filleted section, from linear taper to true flat
  ftaper = f - Rc * (1. - math.cos(theta))
      # depth at which linear taper joins fillet
  LtaperLin = (ftaper/f) * Ltaper
      # length of linear taper that remains 
  if (Lflat == 0.) :
    Ltot = LflatTrue + 2.*(LtaperLin + Lfillet)
    Lflat = Ltot - 2.*Ltaper
      # nominal length of flat if fillet radius Rc were 0
  elif (LflatTrue == 0.) :
    Ltot = Lflat + 2.*Ltaper
    LflatTrue = Ltot - 2.*(LtaperLin + Lfillet)
      # length of the true flat after allowing for fillet
  # print "LtaperLin, Lfillet, LflatTrue = %.4f, %.4f, %.4f" % (LtaperLin,Lfillet,LflatTrue)
  return { "r" : r,  "f" : f, "ang" : ang, "Lflat" : Lflat, "LflatTrue" : LflatTrue, \
     "Ltaper" : Ltaper, "LtaperLin" : LtaperLin, "ftaper" : ftaper, "Rc" : Rc}
 

# -----------------------------------------------------------------------------------

def neatList( pol ) :
  for p in pol :
    print "\n r: %.5f" % p["r"]
    print " f: %.5f" % p["f"]
    print " ang: %.3f" % p["ang"]
    print " Lflat: %.5f" % p["Lflat"]
    print " LflatTrue: %.5f" % p["LflatTrue"]
    print " Ltaper: %.5f" % p["Ltaper"]
    print " LtaperLin: %.5f" % p["LtaperLin"]
    print " ftaper: %.5f" % p["ftaper"]
    print " Rc: %.5f" % p["Rc"]

# -----------------------------------------------------------------------------------

# compute Lflat for retarder that will have phase shift xx at freq ..

def calcLflat() :
  noFlat = oneSection( ang=0., LflatTrue=0., r=r, f=f, Ltaper=Ltaper, Rc=Rc )
    # retarder with no flat at all, just the taper and fillet sections
  vec = propagate( noFlat, Y, "fromOMT" )
  ampPhs(vec)
    # compute amplitude and phase of signal

# -----------------------------------------------------------------------------------

# compute leakage from 70-120 GHz

def computePLeak( pol, outfile="PLeak.dat" ) :
  ofile = open( outfile, "w")
  for fGHz in arange(70.,121.,0.1) :
    v1 = Y 
    if len(pol) > 0 :
      for retarder in pol :
        v2 = propagate( retarder, v1, fGHz ) 
        v1 = v2
    [Yamp,Yphs,Xamp,Xphs,phsdif] = ampPhs(v1)
    Rleak = abs(dot(v1, R)) 		# amplitude of single complex number 
    Lleak = abs(dot(v1, L))
    ofile.write("%6.1f %8.5f %8.5f %8.3f %8.3f %8.3f %8.3f %8.3f\n" % \
       (fGHz, Rleak, Lleak, Yamp, Yphs, Xamp, Xphs, phsdif ) )
  ofile.close()


# -----------------------------------------------------------------------------------
# creates (original) nominal 3 section design

def pnom() :
  pol = []
  pol.append( oneSection( ang=7.50,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  pol.append( oneSection( ang=34.30, LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  pol.append( oneSection( ang=98.55, LflatTrue=0.1250, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  return pol

def newnom() :
  pol = []
  pol.append( oneSection( ang=7.50,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  pol.append( oneSection( ang=34.30, LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  pol.append( oneSection( ang=98.55, LflatTrue=0.1220, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  return pol

# -----------------------------------------------------------------------------------

# propagate vector field through a single retarder section
# inputs:
#   pol = dictionary describing this retarder
#   v1 = vector describing input field
#   direction = 'fromOMT' or 'fromHorn'
# output: vector describing output field

def propagate( pol, v1, fGHz, direction='fromOMT' ) :
  angle = pol["ang"]
  if direction == 'fromHorn' :
    angle = -1. * angle              # angle is reversed if propagating in reverse direction
  v2 = Jrot( v1, angle )             # rotate into reference frame of the retarder
  phfillet = dphiFillet ( r=pol["r"], fstart=pol["f"], fstop=pol["ftaper"], Rc=pol["Rc"], fGHz=fGHz )
  phtaper = dphiTaper( r=pol["r"], f=pol["ftaper"], Ltaper=pol["LtaperLin"], fGHz=fGHz ) 
  phshift = 2. * ( phfillet + phtaper ) + dphi( r=pol["r"], f=pol["f"], L=pol["LflatTrue"], fGHz=fGHz ) 
  v3 = Jdelay( v2, phshift )        # apply phase shift
  return Jrot( v3, -1.*angle )      # rotate back to original frame 

# -----------------------------------------------------------------------------------
             
# return dimension = (xtarg +/- random error)
# error is gaussian distributed with sigma=xtol, but cannot exceed xtol

def dimension( xtarg, xtol ) :
  if (xtol == 0.) : return xtarg
  x = xtarg + 2.*xtol                  # enter loop with dummy value guaranteed to be unacceptable
  while (abs(x - xtarg) > xtol) :
    x = random.gauss( xtarg, xtol )
  return x

# -----------------------------------------------------------------------------------

# add Gaussian errors to polarizer dimensions
# read nominal polarizer list (of retarder dictionaries), add Gaussian errors, return new pol

def addGaussianError( polin, rtol, ftol, angtol, Ltol ) :
  pout=[]
  rErr = dimension( 0., rtol )		# r tolerance applies to all retarder sections  
  if len(polin) == 0 :
    print "error: input polarizer list contains no retarder sections"
  else :
    for rtdr in polin :
      pout.append( oneSection( ang=dimension( rtdr["ang"], angtol ), \
        LflatTrue=dimension( rtdr["LflatTrue"], Ltol ), r=rtdr["r"] + rErr, \
        f=dimension( rtdr["f"], ftol ), Ltaper=rtdr["Ltaper"], Rc=rtdr["Rc"] ) )
  return pout 

# -----------------------------------------------------------------------------------

# leakage of idealized polarizers

def leakIdeal() :

  polA = []				# simple 1 section polarizer
  polA.append( oneSection( ang=45.,  LflatTrue=0.1250, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  computePLeak( polA, outfile="PLeakA.dat" ) 

  polB = []			    # 2 section polarizer, maximally flat, from Kovac
  polB.append( oneSection( ang=15.,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  polB.append( oneSection( ang=75.,  LflatTrue=0.1250, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  computePLeak( polB, outfile="PLeakB.dat" ) 

  polC = []             # 2 section polarizer, CARMA 1mm design
  polC.append( oneSection( ang=15.,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  polC.append( oneSection( ang=74.5,  LflatTrue=0.1250, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  computePLeak( polC, outfile="PLeakC.dat" ) 

  polD = []             # 3 section polarizer, maximally flat, from Kovac
  polD.append( oneSection( ang=6.05,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  polD.append( oneSection( ang=32.35,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  polD.append( oneSection( ang=99.94,  LflatTrue=0.1250, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  computePLeak( polD, outfile="PLeakD.dat" ) 

  polE = []             # 3 section polarizer from Pancharatnam
  polE.append( oneSection( ang=6.50,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  polE.append( oneSection( ang=34.12,  LflatTrue=0.3170, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  polE.append( oneSection( ang=100.69,  LflatTrue=0.1250, f=0.015, Ltaper=0.08, Rc=0.07 ) ) 
  computePLeak( polE, outfile="PLeakE.dat" ) 

  polF = pnom()         # nominal CARMA 3mm design found with blindSearch2
  computePLeak( polF, outfile="PLeakF.dat" ) 

  polG = pnom()	        # same, but with slightly shorter 90 degree retarder
  polG[2]["LflatTrue"] = polF[2]["LflatTrue"] - .003
  computePLeak( polG, outfile="PLeakG.dat" ) 
  
  polH = pnom()	        # same, but with slightly longer 90 degree retarder
  polH[2]["LflatTrue"] = polF[2]["LflatTrue"] + .003
  computePLeak( polH, outfile="PLeakH.dat" ) 

# -----------------------------------------------------------------------------------

# start with nominal design, add Gaussian errors, write out array of leakages vs freq, also avg

def computeLeakage( polin, ntrials=1, rtol=0.0003, ftol=0.0003, angtol=0.1, Ltol=.001 ) :

  freq = arange(70.,121.,1.)
  leakage = zeros( [len(freq),ntrials], dtype=float )    # create array to hold leakages
  effic = zeros( [len(freq),ntrials], dtype=float )      # create array to hold correlation efficiencies
  vsave = []
  pol = []                                               # make pol a global variable!

  for n in range(0,ntrials) :
    if n == 0 :
      pol = polin                                        # first trial is always nominal design, no errors
    else :
      print rtol, ftol, angtol, Ltol
      pol = addGaussianError( polin, rtol, ftol, angtol, Ltol ) 
    m = 0                                                # m is the frequency index
    for fGHz in freq :
      v1 = Y 
      for retarder in pol :
        v2 = propagate( retarder, v1, fGHz ) 
        v1 = v2
      leakage[m,n] = abs(dot(v1, R)) 		          # amplitude of single complex number; note R = L* 
      if (n == 0) :
        vsave.append( array( [v1[0].conjugate(), v1[1].conjugate()] ) )
      effic[m,n] = abs(dot(vsave[m],v1))
      m = m + 1
    print leakage

  ofile1 = open( "leak", "w" )
  ofile2 = open( "effic", "w" )
  m = 0
  for fGHz in freq :
    ofile1.write("%6.1f" % fGHz)
    ofile2.write("%6.1f" % fGHz)
    sum = 0.
    for n in range(ntrials) :
      ofile1.write(" %6.4f" % leakage[m,n])
      sum = sum + leakage[m,n]
      ofile2.write(" %6.4f" % effic[m,n] )
    avg = sum/ntrials
    var = 0.
    for n in range(ntrials) :
      dif = leakage[m,n] - avg
      var = var + dif*dif
    if (ntrials > 1) :
      var = var/(ntrials-1)
    ofile1.write("   %6.4f %6.4f\n" % (avg,sqrt(var)))
    ofile2.write("\n")
    m = m + 1
  ofile1.close()     
  ofile2.close()

# -----------------------------------------------------------------------------------

# compares simulated phase shifts of straight sections with HFSS.04jan2013.csv
# the straight sections use linear tapers with fillet radius of 0

def cmpHFSSstraight( ):
  pol = []
  pol.append( oneSection( ang=0., LflatTrue=0.13727, r=0.0610, f=0.015, Ltaper=0.080, Rc=0. ) ) 
  pol.append( oneSection( ang=0., LflatTrue=0.33494, r=0.0610, f=0.015, Ltaper=0.080, Rc=0. ) )
  pol.append( oneSection( ang=0., LflatTrue=0.50000, r=0.0610, f=0.015, Ltaper=0.080, Rc=0. ) ) 
    # create retarder dictionaries for the 3 retarders simulated with HFSS
  fin = open( "HFSS.04jan2013.csv", "r" )
  fout = open( "HFSS.04jan2013.cmp2.dat", "w" )
  for line in fin:
    a = line.split(',')
    fGHz = float( a[0] )
    fout.write("%7.2f" % fGHz)
    for n in range (0,3) :
      phHFSS = float( a[1 + 3*n] )
      phtaper = dphiTaper( r=pol[n]["r"], f=pol[n]["ftaper"], Ltaper=pol[n]["LtaperLin"], fGHz=fGHz ) 
      ph = 2.*phtaper  + dphi( r=pol[n]["r"], f=pol[n]["f"], L=pol[n]["LflatTrue"], fGHz=fGHz ) 
      fout.write("  %8.3f %8.3f %7.3f" % ( phHFSS, ph, phHFSS-ph ) )
    fout.write("\n")
  fin.close()
  fout.close()

# -----------------------------------------------------------------------------------

# read HFSS csv files, extract relevant columns, compute amps, phases, leakages

def extractHFSS2( ) :
  phsY = []		
  phsX = []
  freq = []
  fin = open( "sphase.csv", "r" )
  for line in fin:
    if line.startswith('"Freq') :
      a = line.split('","')
      for n in range (0, len(a)) :
        if a[n] == "ang_deg(S(2:1,1:1)) [deg]" : Ycol = n
        if a[n] == "ang_deg(S(2:2,1:1)) [deg]" : Xcol = n
      print " Ycol = %d, Xcol = %d" % (Ycol,Xcol)
    else :
      a = line.split(',')
      freq.append( float( a[0] ) )
      phsY.append( float( a[Ycol] ) )
      phsX.append( float( a[Xcol] ) )
      print a[0], a[Ycol], a[Xcol]
  fin.close
  dBY = []
  dBX = []
  fin = open( "smag.csv", "r" )
  m = 0
  for line in fin:
    if line.startswith('"Freq') :
      a = line.split('","')
      for n in range (0, len(a)) :
        if a[n] == "dB(S(2:1,1:1)) []" : Ycol = n
        if a[n] == "dB(S(2:2,1:1)) []" : Xcol = n
      print " Ycol = %d, Xcol = %d" % (Ycol,Xcol)
    else :
      a = line.split(',')
      if float( a[0] ) != freq[m] :
        print "error - freq mismatch; %s, %.3f" % ( a[0],freq[m] )
      m = m + 1
      dBY.append( float( a[Ycol] ) )
      dBX.append( float( a[Xcol] ) )
      print a[0], a[Ycol], a[Xcol]
  fin.close()
  fout = open( "HFSS.dat", "w" )
  for m in range(0, len(freq)) :
    ampY = math.sqrt(pow(10.,dBY[m]/10.))		# convert amps to VOLTAGE
    ampX = math.sqrt(pow(10.,dBX[m]/10.))		
    phdifDeg = phsY[m]-phsX[m]
    if phdifDeg > 180 :  phdifDeg = phdifDeg - 360.
    phdifRad = math.pi * phdifDeg/180.          # convert phase dif to RADIANS
    vec = array( [ (ampX + 1j*0), 
       (ampY*math.cos(phdifRad) + 1j*ampY*math.sin(phdifRad)) ] )
    Rleak = abs(dot(vec, R)) 		             # amplitude of single complex number 
    Lleak = abs(dot(vec, L))
    fout.write("%8.3f   %7.3f %8.3f    %7.3f %8.3f    %8.5f %8.5f %8.3f    %6.4f %6.4f\n" % \
        (freq[m],dBY[m],phsY[m],dBX[m],phsX[m],ampY,ampX,phdifDeg,Rleak,Lleak ) )
  fout.close()

# absolute phase shift through 0.7" long retarder + 0.3" long regular guide
# meant to compare with simulation of SlowWaveRetarder

def cmpSlow() :
  fout = open( "cmpSlow2.dat", "w" )
  L0 = 0.85   # length of circular guide
  L1 = 0.15   # length of faceted guide
  [ fcX0, fcY0 ] = cutoff( .061, .0 )    
  [ fcX1, fcY1 ] = cutoff( .061, .015/.061 )    
  for fGHz in arange( 70., 120.1, .1 ) :
    phiX = -360. * 2.54 / clight * ( L0 * math.sqrt( pow(fGHz,2) - pow(fcX0,2) ) \
                                  + L1 * math.sqrt( pow(fGHz,2) - pow(fcX1,2) ) )
    phiY = -360. * 2.54 / clight * ( L0 * math.sqrt( pow(fGHz,2) - pow(fcY0,2) ) \
                                  + L1 * math.sqrt( pow(fGHz,2) - pow(fcY1,2) ) )
    phdif = phiX - phiY
    phiX = phiX % 360.
    if phiX > 180. : phiX = phiX - 360.
    phiY = phiY % 360.
    if phiY > 180. : phiY = phiY - 360.
    fout.write("%8.2f  %10.3f  %10.3f  %10.3f\n" % ( fGHz, phiX, phiY, phdif ) )
  fout.close()

# extracts differential phase shift through RETARDER section, and reflection coefficients,
#   from HFSS csv files
def extractHFSS( HFSSphs, HFSSmag, outfile ) :
  phsY = []		
  phsX = []
  freq = []
  fin = open( HFSSphs, "r" )
  print HFSSphs
  for line in fin:
    if line.startswith('"Freq') :
      a = line.split('","')
      for n in range (0, len(a)) :
        if a[n] == "ang_deg(S(2:1,1:1)) [deg]" : Ycol = n
        if a[n] == "ang_deg(S(2:2,1:2)) [deg]" : Xcol = n
      print " Ycol = %d, Xcol = %d" % (Ycol,Xcol)
    else :
      a = line.split(',')
      if len(a[1]) > 0 :       # avoids lines at end > 120 GHz with no data
        freq.append( float( a[0] ) )
        phsY.append( float( a[Ycol] ) )
        phsX.append( float( a[Xcol] ) )
        # print a[0], a[Ycol], a[Xcol], float(a[Ycol])-float(a[Xcol])
  fin.close()

  S11Y = []
  S11X = []
  fin = open( HFSSmag, "r" )
  print HFSSmag
  m = 0
  for line in fin:
    if line.startswith('"Freq') :
      a = line.split('","')
      for n in range (0, len(a)) :
        if a[n] == "dB(S(1:1,1:1)) []" : Ycol = n
        if a[n] == "dB(S(1:2,1:2)) []" : Xcol = n
      print " Ycol = %d, Xcol = %d" % (Ycol,Xcol)
    else :
      a = line.split(',')
      if len(a[1]) > 0 :       # avoids lines at end > 120 GHz with no data
        if float( a[0] ) != freq[m] :
          print "error - freq mismatch; %s, %.3f" % ( a[0],freq[m] )
        m = m + 1
        S11Y.append( float( a[Ycol] ) )
        S11X.append( float( a[Xcol] ) )
      #print a[0], a[Ycol], a[Xcol]
  fin.close()

  fout = open( outfile, "w" )
  for m in range(0, len(freq)) :
    phdif = phsX[m]-phsY[m]
    if phdif < 0. : phdif = phdif + 360.
    fout.write( "%8.2f  %11.7f  %11.7f  %11.7f  %10.3f  %10.3f\n" % \
      (freq[m],phsY[m],phsX[m],phdif,S11Y[m],S11X[m]) )
  fout.close()
  
# for slow wave leakage calculations, read ret90 and ret180 phase shifts from table
# ang180a = angle of 1st 180 degree retarder section (skip if 0)
# ang180b = angle of 2nd 180 degree retarder section (skip if 0)
# ang90 = angle of 90 degree retarder section

def computeSWLeak( ang180a, ang180b, ang90, outfile ) :
  ntrials = 5
  freq = []
  ret90 = []
  ret180 = []
  offset= [-2.,-1.5,-1.,-0.5,0.,0.5,1.,1.5,2.]

  fin = open( "SW3+6.dat", "r" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      freq.append( float(a[0]) )
      #ret90.append( float(a[9]) )
      ret180.append( float(a[3]) )
      ret90.append( ret180[ len(ret180)-1 ]/2. )
      print freq[len(freq)-1], ret90[len(ret90)-1], ret180[(len(ret180)-1)]
  fin.close()

  leakage = zeros( [len(freq),ntrials], dtype=float )    # create array to hold leakages
  effic = zeros( [len(freq),ntrials], dtype=float )      # create array to hold correlation efficiencies
  ofile1 = open( outfile, "w" )
  

  m = 0                                                # m is the frequency index
  for fGHz in freq :
    ofile1.write("%6.1f" % fGHz)
    for off in offset:
      v1 = Y 
      if ang180a != 0. :
        v2 = Jrot( v1, ang180a )            # rotate into reference frame of the retarder
        v3 = Jdelay( v2, ret180[m] )        # apply phase shift
        v1 = Jrot( v3, -1.*ang180a )        # rotate back to original frame 
      if ang180b != 0. :
        v2 = Jrot( v1, ang180b )            # rotate into reference frame of the retarder
        v3 = Jdelay( v2, ret180[m] )        # apply phase shift
        v1 = Jrot( v3, -1.*ang180b )        # rotate back to original frame 
      if ang90 != 0. :
        v2 = Jrot( v1, ang90+off )            # rotate into reference frame of the retarder
        v3 = Jdelay( v2, ret90[m] )        # apply phase shift
        v1 = Jrot( v3, -1.*(ang90+off) )        # rotate back to original frame 
      leakage = abs(dot(v1, R)) 		          # amplitude of single complex number; note R = L* 
      ofile1.write(" %6.4f" % leakage)
    ofile1.write("\n")
    m = m + 1
  ofile1.close()

def doit() :
  extractHFSS( "SW3phase.csv", "SW3mag.csv", "SW3.dat" )
  extractHFSS( "SW6_phase.csv", "SW6_mag.csv", "SW6.dat" )
