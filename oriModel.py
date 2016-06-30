# model.py
# misc python calculations relevant to Orion HII region paper

# from Numeric import *
from numpy import *
import math

# ------------- constants ------------------ #
Te = 8000            # electron temp in K
kB = 1.3806e-16       # boltzmann's constant, erg K-1
h = 6.626e-27        # plank's constant, erg s
c = 2.9979e10        # c in cm/sec
clight = 2.9979e10   # c in cm/sec
au = 1.496e13        # 1 AU in cm
pc = 3.0856e18		 # 1 pc in cm
G = 6.6732e-8        # grav constant in dynes cm^2 g^-2
dAU = 415 * pc/au    # 400 pc in AU 
mass_e = 9.109e-28	 # electron mass in grams
mass_H = 1.6605e-24  # atomic mass unit in grams
mass_sun = 1.989e33  # solar mass in g
SahaFactor = pow( 2. * math.pi * mass_e * kB ,1.5) / pow( h, 3. )

def boxcargauss( v1, vstep, nchan, a0, a1, v0, fwhm, outfile ) :
  ofile = open(outfile,"w")
  vcenter = v1
  amp = 0.
  nsteps = 0

  # --- step through velocity range in 0.1 km/sec steps --- #
  for v in arange( v1-vstep/2., v1+nchan*vstep, .1 ) :
    if (v > (vcenter+vstep/2.) ) :
      ofile.write("%10.2f  %10.4f  %d\n" % (vcenter,amp/nsteps,nsteps))
      amp = 0.
      nsteps = 0
      vcenter = vcenter + vstep
    argument = (v-v0)*(v-v0)/(fwhm*fwhm)
    if (argument > 100.) :
      amp = amp + a0
    else :
      amp = amp + a0 + a1 * exp(-2.772 * argument)
    nsteps = nsteps + 1
  ofile.close()

def halpha() :
#  boxcargauss(-80.,4.,44,.09,.047,23.8,22.6,"h41mod4")
#  boxcargauss(-300.,50,20,.09,.047,23.8,22.6,"h41mod50")
#  boxcargauss(-681.8,40.4,35,.175,.55,23.8,22.6,"h30mod40")
#  boxcargauss(-681.8,1,1400,.175,.55,23.8,22.6,"h30mod1")

#  boxcargauss( v1, vstep, nchan, a0, a1, v0, fwhm, outfile ) 
#  boxcargauss(  -60, 1, 151, .0286, .0104, 20.1, 39.0, "h53mod1" ) 

#  boxcargauss( -200, 100, 5, 0.175, 0.550, 0., 22.6, "h30flux")
#  boxcargauss( -200, 100, 5, 0.090, 0.047, 0., 22.6, "h41flux")
#  boxcargauss( -200, 100, 5, .0286, .0104, 0., 39.0, "h53flux")

#   boxcargauss( -200., 100., 5, .25, .15, 23.8, 25, "h30Imod100" )  # used to check intensity in 100 km/sec chan
  boxcargauss( -1005,  40.4, 43, .175, .55, 23.8, 22.6,"h30mod40")
  boxcargauss( -1005., 40.4, 43, .250, .15,  5,   25,  "h30Imod40" ) 
  
  
# --- compute flux density of optically thick region of brightness temp Tb, radius r AU, temp Tb
def spectrum(fGHz, Tb, rAU=10 ) :
  # begin by computing the reference flux (optically thick, radius 10 AU, temp T)
  d = 415 * 3.0856e18  # 415 pc in cm
  r0 = rAU * 1.496e13      # 10 AU in cm
  #Te = 8200               # electron temp in K
  kB = 1.3806e-16          # boltzmann's constant
  c = 2.9979e10           # c in cm/sec
  S0 = math.pi * r0*r0 * 2.*kB*Tb * (fGHz * 1.e9) * (fGHz * 1.e9) / (d * d * c * c * 1.e-23)
  print "source radius = %.2f AU at 415 pc" % rAU
  print "brightness temp = %.0f K" % Tb
  print "-> flux density = %.4f mJy" % (S0*1.e3)
  return S0

  #tau0 = 8.235e-2 * pow(Te,-1.35) * pow(fGHz,-2.1)	# Lang, eq 1-223
  
# --- continuum and line emissivities vs radius --- #
def radial (outfile, fGHz, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) :
  ofile = open(outfile,"w")
  ofile.write("# alpha = %.2f\n" % alpha)
  ofile.write("# r0 = %.2f AU\n" % r0AU )
  ofile.write("# n0 = %.2e cm-3 at r0\n" % n0)
  ofile.write("# rmax = %.2e AU\n" % rmaxAU)
  ofile.write("# dx = %.4e AU\n" % dxAU)
  ofile.write("# f = %.2f GHz\n" % fGHz)
  for xAU in arange(dxAU/2.,rmaxAU,dxAU) :
    EM = em(xAU, alpha, r0AU, n0, rmaxAU)
    tC = tauC(EM,fGHz)
    tL = tauL(EM,fGHz,fwhm)
    eC = 1.
    if (tC < 100.): 
      eC = (1. - exp(-1.*tC))
    eL = 1.
    if ((tL+tC) < 100.): 
      eL = 1. - exp(-1.*(tL+tC))
    ofile.write("%8.3f  %8.5f  %8.5f\n" % (xAU, eC, eL))
  ofile.close()

# --- continuum opacity, Rohlfs & Wilson 9.35 --- #
def tauC( EM, fGHz ) :
  return (8.235e-2 * pow(Te, -1.35) * pow(fGHz, -2.1) * EM)

# --- line opacity, Rohlfs & Wilson 13.27 --- #
def tauL( EM, fGHz, fwhm ) :
  dfkHz = 1.e11 * fwhm * fGHz / c 
  return (1.92e3 * pow(Te, -2.5) * EM / dfkHz )

# --- find emission measure at projected distance x(cm) from the center --- #
def em (xAU, alpha, r0AU, n0, rmaxAU) :
  dyAU = .1
  ymaxAU = math.sqrt(rmaxAU*rmaxAU - xAU*xAU)		# integrate to same max radius each time
  em = 0.
  for yAU in arange( dyAU/2., ymaxAU, dyAU ) :
    rAU = math.sqrt(xAU * xAU + yAU * yAU)
    ne = n0
    if (rAU > r0AU) :
      ne = n0 * pow(rAU/r0AU, alpha)	# electron density at this radius, cm-3
    em = em + ne * ne * dyAU * au/pc    # emission measure in cm-6 pc
  return 2.*em							# x 2 because we integrated through half the column  

# --- compute flux densities in MILLIJANSKY for continuum and at peak of line vs freq --- #
def flux( outfile, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) :
  flist = [43.,90.,231.9, 353.6, 662.4 ]
  ofile = open(outfile, "w" )
  ofile.write("# r0 = %.2f AU\n" % r0AU)
  ofile.write("# n0 = %.2e cm^-3\n" % n0)
  ofile.write("# alpha = %.2f (ne index)\n" % alpha)
  ofile.write("# rmax = %.2f AU\n" % rmaxAU)
  ofile.write("# dx = %.3f AU\n" % dxAU)
  ofile.write("# fwhm = %.2f km/sec\n" % fwhm)
  const = 4. * math.pi * kB * Te * 1.e44 / (dAU * dAU * c * c)
  #for fGHz in flist :
  for logf in arange(0.4, 3.6, 0.1):
    fGHz = pow(10., logf)
    print fGHz
    logf = math.log10(fGHz)
    sumC = 0.
    sumL = 0.
    for xAU in arange(dxAU/2.,rmaxAU,dxAU) :
      EM = em(xAU, alpha, r0AU, n0, rmaxAU)    # emission measure in annulus at radius xAU

      tC = tauC(EM,fGHz)					   # continuum opacity for this emission measure
      eC = 1.
      if (tC < 100.): 
        eC = (1. - exp(-1.*tC))	               # 1-exp(-tauC) for this annulus
      sumC = sumC + eC * xAU * dxAU            # add (1-exp(-tauC)) * radius * dr

      tL = tauL(EM,fGHz,fwhm)		           # opacity at center of gaussian line
      eL = 1.
      if ((tL+tC) < 100.):
        eL = (1. - exp(-1.*(tL+tC)))
      sumL = sumL + eL * xAU * dxAU
      
      # --- the following section finds average opacity of line in wider channel --- #
      # eL = 0.
      # ntmp = 0
      # for v in arange(2.5,47.5,5.) :	# integrate from 0-50 km/sec in 5 km/sec bins
      #  tLv = tL * gauss(v,fwhm)
      #  eLv = 1.						# eLv is (1-exp(-tau)) at this velocity
      #  if ((tLv + tC) < 30) :
      #    eLv = 1. - exp(-1.*(tLv + tC))
      #  eL = eL + eLv
      #  ntmp = ntmp + 1
      # eL = eL/ntmp						# eL is average eLv within the 100 km/sec wide chan
      # sumL = sumL + eL * xAU * dxAU

    SC = const * sumC * fGHz*fGHz
    SL = const * sumL * fGHz*fGHz       
    ofile.write("%8.3f  %8.4f  %8.3f  %8.4f  %8.3f  %8.4f\n" %  \
       (fGHz,logf,SC,math.log10(SC),SL,math.log10(SL)))

# --- find gauss at velocity v relative to velocity 0, for width fwhm --- #
def gauss( v, fwhm ) :
  return exp(-2.7726 * pow(v/fwhm, 2.))

def examples () :
# radial (outfile, fGHz, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) 
#  radial ("r43.dat", 43., 39., -10., 5., 1.e7, 50., 0.5 ) 
# flux ( outfile,     fwhm, alpha,  r0AU,   n0,  rmaxAU,  dxAU ) :
# flux ("fnormal.dat", 26.,  -10.,  10.0,  2.e7,   20.,   0.5 )
#  flux ("BNmodel.dat",      23.,  -5.,    9.5,  3.e7,   30.,   0.5 )   # save
#  flux ("BNmodel1.dat",      23.,  -3.3,    7,  4.5e7,   30.,   0.5 )   # good fit to most of the spectrum
#  flux ("HII.dat",     23.,  -20.,   9.5,  3.e7,   30.,   0.5 )
# def flux (outfile, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) :
#  flux ("Imodel.dat",   23.,  -10,   6.8,  8.5e7,   60.,   0.5 )   # save - good fit
#  flux ("IHII.dat",   23.,  -100,   7.6,  7.1e7,   60.,   0.5 )   # save - good fit
#  flux ("BNmodel2.dat",      23.,  -10.,    10.5,  2.8e7,   30.,   0.5 )   # good fit to 43,89,229 GHz points
#  flux ("Imodel.dat",   23.,  -10,   6.8,  8.5e7,   60.,   0.5 )   # save - good fit
#  flux ("classicHII.dat",   23.,  -20,   6.8,  1.e7,   60.,   0.5 ) 
#  flux ("classicWind.dat",   23.,  -2,   .1,  1.e10,   60.,   0.5 )
#  flux ("classicHyper.dat",   23.,  -3,   3,  1.e7,   60.,   0.5 )

#   ---- feb 2011 ---
#   flux ("Imodel.26feb2011.dat",   23.,  -20, 6.8,  3e8,  50., 0.5 )

#   ---- jul 2012 ----
#   radial(outfile,           fGHz, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) 
#   radial("Imodel.radial.dat", 662., 23., -50., 6.8,   3e8,  20,  0.1)
#   flux ("Imodel.dat",   23.,  -20, 6.8,  3e8,  10., 0.1 )
#   flux ("Imodelextreme.dat", 40., -50., 6.8, 6e8, 10., 0.1 )
#   flux ("Imodel3e8.dat",   23.,  -50, 6.8,  3e8,  10., 0.1 )
#   flux ("Imodel3e9.dat",   23.,  -50, 6.8,  3e9,  10., 0.1 )
#   flux ("Imodel3e10.dat",   23.,  -50, 6.8,  3e10,  10., 0.1 )

#   flux ( outfile,     fwhm, alpha,  r0AU,   n0,  rmaxAU,  dxAU ) :
   flux ("Imodel.dat",   30.,  -50,   7.5,  1.5e8,   20.,   0.1 )   # save - good fit
   flux( "BNmodel2.dat",  30.,  -3.5,  7.4,  5e7,   50.,   0.1 )	# good fit feb2011

# freq for Halpha(n+1 -> n); for H and heavy element (K); velocity difference between them
def freq( n ):
  RH = 3.28805129e6	 # GHz
  RK = 3.28984202e6	 # GHz; for heavy atom like K
 
  fK = (RK * (1./float(n*n) - 1./float((n+1)*(n+1))))
  print fH, fK, -3.e5*(fK-fH)/fH

# --- eqn 20 of Dalgarno and Lane --- #
def I1( a ) :
  x = .1 
  y = 1000.
  while (y > 0.05) :
    x = 2. * x
    y = exp(-1. * a * x) * (pow(x,0.5)*pow(x+1.,1.5) + pow(x,1.5)*pow(x+1.,0.5))
  xmax = x		# this value of x is after function y(x) tapers back toward zero
  dx = xmax/1000.     # 1000 steps
  tot = 0.
  for x in arange(dx/2., xmax, dx) :
    y = exp(-1. * a * x) * (pow(x,0.5)*pow(x+1.,1.5) + pow(x,1.5)*pow(x+1.,0.5))
    tot = tot + y * dx
  tot = tot * (1. - exp(-1.*a))
  return tot

def I2( a ) :
  x = .1 
  y = 1000.
  while (y > 0.05) :
    x = 2. * x
    y = exp(-1.*a*x) * (x*pow(x+1.,1.5) + pow(x,1.5)*(x+1.))
  xmax = x		# this value of x is after function y(x) tapers back toward zero
  dx = xmax/1000.     # 1000 steps
  tot = 0.
  for x in arange(dx/2., xmax, dx) :
    y = exp(-1.*a*x) * (x*pow(x+1.,1.5) + pow(x,1.5)*(x+1.))
    tot = tot + y * dx
  tot = tot * (1. - exp(-1.*a))
  return tot

def I3( a ) :
  x = 1. 
  y = 1000.
  while (y > 0.05) :
    x = 2. * x
    y = exp(-1.*a*x) * (math.log(x) + math.log(x+1.)) * pow(x,1.5) * pow(x+1.,1.5)
  xmax = x		# this value of x is after function y(x) tapers back toward zero
  dx = xmax/1000.     # 1000 steps
  tot = 0.
  for x in arange(dx/2., xmax, dx) :
    y = exp(-1.*a*x) * (math.log(x) + math.log(x+1.)) * pow(x,1.5) * pow(x+1.,1.5)
    tot = tot + y * dx
  tot = tot * (1. - exp(-1.*a))
  return tot

def I4( a ) :
  x = .1 
  y = 1000.
  while (y > 0.05) :
    x = 2. * x
    y = exp(-1.*a*x) * (pow(x,2.)*pow(x+1.,1.5) + pow(x,1.5)*pow(x+1.,2.))
  xmax = x		# this value of x is after function y(x) tapers back toward zero
  dx = xmax/1000.     # 1000 steps
  tot = 0.
  for x in arange(dx/2., xmax, dx) :
    y = exp(-1.*a*x) * (pow(x,2.)*pow(x+1.,1.5) + pow(x,1.5)*pow(x+1.,2.))
    tot = tot + y * dx
  tot = tot * (1. - exp(-1.*a))
  return tot

def kv( D1, D2, D3, D4, D5, T, fGHz) :
  e = 4.1344e-6 * fGHz   # 12395/lambda_Angstroms
  a = 4.8014e-2 * fGHz/T    # hv/kT
  term = D1 * I1(a) + pow(e,0.5)* D2 * I2(a) + e * D3 * I3(a) \
         + e * D3 * 3/a * I1(a) * math.log(e) + e * D4 * 3./a * I1(a) \
         + pow(e,1.5) * D5 * I4(a)
  return 3.269e-18 * pow(T,-2.5) * term

def kHI ( T, fGHz ) :
  D1 = 11.240
  D2 = 7.200
  D3 = 4.959
  D4 = -11.580
  D5 = -0.4023
  return kv( D1, D2, D3, D4, D5, T, fGHz )

def kHIRM ( T ) :
  a0 = 3.376e-17
  a1 = -2.149e-20
  a2 = 6.646e-24
  a3 = -7.853e-28
  return a0 + a1*T + a2*pow(T,2.) + a3*pow(T,3.)

# --- compare my evaluation of eqn 20 of Dalgarno and Lane with fit from Reid and Menten --- #
# --- this shows that the fit is excellent! --- #
def HIcmp () :
  for T in arange(1000.,3000.,50.) :
    print T,1.e18*kHIRM(T), 1.e18*kHI(T,10.)

def kH2 ( T, fGHz ) :
  D1 = 1.6456
  D2 = 5.8492
  D3 = 0.86147
  D4 = -2.6502
  D5 = 0.
  return kv( D1, D2, D3, D4, D5, T, fGHz )

def kH2RM ( T ) :
  a0 = 8.939e-18
  a1 = -3.555e-21
  a2 = 1.100e-24
  a3 = -1.319e-28
  return a0 + a1*T + a2*pow(T,2.) + a3*pow(T,3.)

def H2cmp () :
  for T in arange(1000.,3000.,50.) :
    print T,1.e18*kH2RM(T), 1.e18*kH2(T,10.)


# --- fractional ionization of K at temperature T
def nK ( nH, T, ne ) :
  fK = 1.32e-7
  I = 4.34  # eV  
  C = 2. * SahaFactor * pow(T,1.5) * exp(-I*11606./T) / (ne)
  return fK * nH * C/(1.+C)

# --- fractional ionization of Na at temperature T
def nNa ( nH, T, ne ) :
  fNa = 1.887e-6
  I = 5.14  # eV  
  C = 2. * SahaFactor * pow(T,1.5) * exp(-I*11606./T) / (ne)
  return fNa * nH * C/(1.+C)

# --- fractional ionization of Ca at temperature T
def nCa ( nH, T, ne ) :
  fCa = 2.267e-6
  I = 6.11  # eV  
  C = 2. * SahaFactor * pow(T,1.5) * exp(-I*11606./T) / (ne)
  return fCa * nH * C/(1.+C)

# --- fractional ionization of Al at temperature T
def nAl ( nH, T, ne ) :
  fAl = 2.67e-6
  I = 5.99  # eV  
  C = 2. * SahaFactor * pow(T,1.5) * exp(-I*11606./T) / (ne)
  return fAl * nH * C/(1.+C)

# --- solve iteratively for electron density at density nH, temperature T --- #
def neSolve ( nH, T ) :
  niter = 1
  neGuess = 1.e-7 * nH
  nIon = nK(nH,T,neGuess ) + nNa(nH,T,neGuess) + nCa(nH,T,neGuess) + nAl(nH,T,neGuess)
  error = abs(neGuess - nIon)/neGuess
  while ((error > .005) & (niter < 20)) :
    neGuess = (neGuess + nIon)/2.
    #print neGuess, nIon, error
    nIon = nK(nH,T,neGuess ) + nNa(nH,T,neGuess) + nCa(nH,T,neGuess) + nAl(nH,T,neGuess)
    error = (neGuess - nIon)/neGuess
    niter = niter + 1
  print "%.0f  %.3e  %.3e  %.3e  %.3e  %.3e" % ( T,nK(nH,T,neGuess ),nNa(nH,T,neGuess),nCa(nH,T,neGuess),nAl(nH,T,neGuess),neGuess )
  return neGuess  

# try to replicate Reid and Menten fig 5
def RMfig5() :
  nH = 1.e12
  for T in arange(1000,2020,20) :
    neSolve( nH, T )
  
# --- number of ionizing photons from eq 2 of Mezger, Smith, & Churchwell 1974 --- #
def Mezger( fGHz, tK, SJy ) :
  Dkpc = 0.41
     # distance to Orion in kpc
  a = 0.366 * pow(fGHz,0.1) * pow(tK,-.15) * ( math.log(4.995e-2/fGHz) + 1.5*math.log(tK) )
     # from Mezger and Henderson, eqn A.2
  Nu = 4.76e48 * pow(fGHz,0.1) * pow(tK,-.45) * SJy * pow(Dkpc,2.) / a
  return Nu

# --- Gaunt factor, eqn 8.27 in Rohlfs --- #
def gaunt( fGHz, Telec ) :
  return math.log( 4.955e-2/fGHz) + 1.5*math.log(Telec)

# --- free-free opacity from eqn 8.26 in Rohlfs --- #
def tauff( fGHz, EM, Telec) :
  return 3.014e-2 * pow(Telec, -1.5) * pow(fGHz, -2.) * EM * gaunt(fGHz, Telec) 

# --- compute minimum mass of region with opacity tau and flux density S at freq fGHz --- #
def minmass( tau, SJy, Te, fGHz ) :
  EM = tau / tauff( fGHz, 1., Te )
  print "EM = %.2e pc cm-6" % EM
  Tb = (1. - math.exp(-tau)) * Te
  print "Tb = %.0f K" % Tb
  S10AU = spectrum( fGHz, Tb )			# flux density in Jy for optically thick emission 10 AU radius 
  print "S10AU = %.2e Jy" % S10AU
  rAU = 10. * math.sqrt(SJy/S10AU)
  print "r = %.3f AU" % rAU
  ne = math.sqrt( EM / (1.496e13 * rAU / 3.086e18) )   # ne cm-3 assuming thickness = radius
  print "ne = %.2e cm-3" % ne
  mass = math.pi * pow(rAU * 1.496e13, 3.) * ne * 1.66e-24 / 1.99e33
  print "total mass = %.2e Mo" % mass

# --- compute accretion luminosity for Mdot solar masses/yr falling onto star of mass M, radius Rinner --- #
def accLum( Mdot=1.e-3, M=10., rAU=0.1 ) :
  print "M = %.2f Mo, Mdot = %.2e Mo/yr, inner radius = %.2f AU" % (M, Mdot, rAU)
  acclum = G * (M * 1.99e33) * (Mdot * 1.99e33 / 3.154e7) / (rAU * au * 3.826e33)
  print "accretion luminosity = %.2e Lo/yr" % acclum

#   flux ( outfile,     fwhm, alpha,  r0AU,   n0,  rmaxAU,  dxAU ) :
#   flux( "BNmodel2.dat",  30.,  -3.3,  7.3,  5.e7,   30.,   0.2 )	# save
# flux( "BNmodel2.dat",  30.,  -3.5,  7.4,  5e7,   40.,   0.1 )	# good fit feb2011
# flux ("Imodel.dat",   30.,  -50,   7.5,  8.5e7,   20.,   0.1 )   # fits SMA point at 348 GHz
# flux ("Imodel.dat",   30.,  -50,   7.5,  1.5e8,   20.,   0.1 )   # save - good fit

# see if I can reproduce Reid and Menten fig 6
# unfortunately, no: need radio of H/H2 in order to do it right
# at 2000 K most of the hydrogen must be H, so possibly I get about the right answer
def RMfig6( ) :
  #fGHz = 10.        # freq in GHz
  #tslab = 3.e12     # thickness of slab in cm
  #nH = 2.e12        # total H density = nH + 2.*nH2 in cm^-3
  fGHz = 43
  tslab = 3.e13
  nH = 1.e11
  for T in arange(1000.,2100.,100.) :
    ne = neSolve( nH, T )
    tau = nH * ne * kB * T * kHI( T, fGHz) * tslab   # ignoring contribution from H2
    print "%.0f  %.3e  %.3f" % ( T, tau, math.log10(tau) )

# to satisfy the referee for SrcI paper, I need to compute ratio of free-free to H- opacity
# I will compute this as a function of temperature and frequency
#   - to compute free-free opacity, assume gas is fully ionized
#   - to compute H- opacity, assume gas is fully neutral
# temp should be > 2000 K to make sure that most H is atomic
def opacityRatio( T, fGHz ):
  nH = 2.e11
  tslab = 3.e12
  ne = neSolve( nH, T )
  EM = nH * nH * tslab/3.08e18 
  tau_ff = tauff( fGHz, EM, T ) 
  tau_Hminus = nH * ne * kB * T * kHI( T, fGHz ) * tslab
  ratio = (tau_ff/nH) / (tau_Hminus/ne)
  print " tau_ff = %.2e, tau_Hminus = %.2e, ratio per electron = %.2e" % (tau_ff,tau_Hminus,ratio)

# try to see turnover in Hminus spectrum
# compute flux density in milliJy!
# input parameters: T = temperature in K, nH = total hydrogen density nH+2nH2 in cm-3,
#   and tslab = thickness of slab in cm
def HminusSpectrum( T, nH, radiusAU, tslabAU ):
  radius = radiusAU		# radius of slab in AU
  tslab = 1.5e13 * tslabAU
  fout = open( "Hminus.dat", "w" )
  fout.write("# T = %.1f K\n" % T)
  fout.write("# nH = %.3e cm^-3\n" % nH)
  fout.write("# slab thickness = %.2f AU\n" % tslabAU)
  fout.write("# slab radius = %.2f AU\n" % radiusAU )

  # Note: on 6/29/2016 I found the following formula, which seems to be high by factor of 2
  # I have now rewritten the formula to be more transparent below, without using 'const'
  # const = 4. * math.pi * kB * 1.e44 * pow( radiusAU/(dAU*c), 2)

  # Note: on 6/29/2016 I found the following formula, which does not include correction for helium and metals
  # totmass = nH * tslab * pow( radius*1.5e13, 2.) * mass_H / mass_sun     

  totmass = 1.36 * nH * tslab * math.pi * pow( radiusAU*1.5e13, 2.) * mass_H / mass_sun     
  fout.write("# tot mass = %.4f Mo\n" % totmass )
  print "totmass = %.4f Mo" % totmass
  ne = neSolve( nH, T )  
  for logf in arange(0.,4.,.1) :
    fGHz = pow(10.,logf)
    tau_Hminus = nH * ne * kB * T * kHI( T, fGHz ) * tslab
    if (tau_Hminus > 10.) :
      Tb = T
    else :
      Tb = T * (1 - math.exp(-1.*tau_Hminus))
    # SmJy = const * Tb * fGHz*fGHz       
    SmJy = 2. * kB * Tb * pow(1.e9*fGHz/clight, 2) * math.pi * pow(radiusAU/dAU, 2.) / 1.e-26
    fout.write("%.1f  %.1f  %.4f  %.4f\n" % (fGHz, Tb, tau_Hminus, SmJy) )
  fout.close()

# same thing, but for regular brehmsstrahlung
# assume slab is fully ionized
def BrehmSpectrum( T, ne, tslab=3.e12 ) :
  radius = 5
  const = 4. * math.pi * kB * 1.e44 * pow( radius/(dAU*c), 2)
  fout = open( "Brehm.dat", "w" )
  fout.write("# T = %.1f K\n" % T)
  fout.write("# ne = %.3e cm^-3\n" % ne)
  fout.write("# slab thickness = %.3e cm\n" % tslab)
  fout.write("# slab radius = %.2f AU\n" % radius )
  totmass = ne * tslab * pow( radius*1.5e13, 2.) * mass_H / mass_sun
  fout.write("# tot mass = %.4f Mo\n" % totmass )
  EM = tslab * ne * ne / pc
  fout.write("# emission meas = %.3e cm^6 pc^-1\n" % EM )
  for fGHz in arange( 10., 700., 10.) :
    tau_brehm = tauff( fGHz, EM, T )
    if (tau_brehm > 10.) :
      Tb = T
    else :
      Tb = T * (1 - math.exp(-1.*tau_brehm))
    SmJy = const * Tb * fGHz*fGHz       
    fout.write("%.1f  %.1f  %.4f  %.4f\n" % (fGHz, Tb, tau_brehm, SmJy) )
  fout.close()
 
   
# compute density of H2 at temperature T, given total H density nH = n(H) + 2*n(H2)     
# use equilibrium constant Kp from Tsuji 1973, table 2
def nH( T, NH ) :
  theta = 5040./T
  logKp = 1.2739e1 - 5.1172e0 * theta + 1.2572e-1 * pow(theta,2.) - 1.4149e-2 * pow(theta,3.) + 6.0321e-4 * pow(theta,4)
  Kp = pow(10.,logKp)
  print Kp
  a = kB * T
  b = Kp/2.
  c = -Kp * NH / 2.
  if ( b*b - 4*a*c) < 0. :
    print "error"
  else :
    nH = -b + math.sqrt( b*b - 4. * a * c )/(2.*a)
    print "nH = %.3e, nH/nH = %.3e" % (nH, (nH/NH) )
  
