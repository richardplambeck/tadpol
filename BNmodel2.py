# BNmodel.py
# adapted from oriModel.py, 2020-05-22
# misc python calculations relevant to Orion HII region paper

import numpy
import math
import pickle
import string
import matplotlib
# matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
from scipy.ndimage import rotate
from scipy.special import kn
from scipy.signal import convolve2d

# ------------- constants ------------------ #
Te = 8000            # electron temp in K
kB = 1.3806e-16       # boltzmann's constant, erg K-1
h = 6.626e-27        # plank's constant, erg s
hplanck = 6.626e-27        # plank's constant, erg s
#c = 2.9979e10        # c in cm/sec
clight = 2.9979e10   # c in cm/sec
au = 1.496e13        # 1 AU in cm
pc = 3.0856e18		 # 1 pc in cm
G = 6.6732e-8        # grav constant in dynes cm^2 g^-2
dAU = 400 * pc/au    # 400 pc in AU 
mass_e = 9.109e-28	 # electron mass in grams
mass_H = 1.6605e-24  # atomic mass unit in grams
mass_sun = 1.989e33  # solar mass in g
SahaFactor = pow( 2. * math.pi * mass_e * kB ,1.5) / pow( h, 3. )
e_esu = 4.08e-10
dist_AU = 400.

# compute flux density in milliJy for blackbody of temperature Tb
def blackbody(fGHz, Tb, dmajArcsec, dminArcsec) :
  dmaj = math.pi * dmajArcsec/(60.*60.*180.)    # convert to radians
  dmin = math.pi * dminArcsec/(60.*60.*180.)    # convert to radians
  dOmegaEllipse = math.pi * dmaj/2. * dmin/2.
  dOmegaGaussian = 1.13 * dmaj * dmin	        # from notebook p. 98
  B = 2. * hplanck * pow(fGHz*1.e9,3)/pow(clight,2) / (math.exp(hplanck*fGHz*1.e9/(kB*Tb)) - 1.) 
  S0 = dOmegaEllipse * B/ 1.e-26
  S0RJ = dOmegaEllipse * 2. * kB * Tb * pow(fGHz*1.e9/clight,2)/ 1.e-26
  S1 = dOmegaGaussian * B/ 1.e-26
  S1RJ = dOmegaGaussian * 2. * kB * Tb * pow(fGHz*1.e9/clight,2) /1.e-26
  print "flux density of %.1f K blackbody uniform disk, %.3f x %.3f arcsec ellipse: %.4f mJy (RJ: %.4f mJy)" % (Tb,dmajArcsec,dminArcsec,S0,S0RJ)
  print "flux density of %.1f K blackbody gaussian, %.3f x %.3f FWHM: %.4f mJy (RJ: %.4f mJy)" % (Tb,dmajArcsec,dminArcsec,S1,S1RJ)
  
# --- compute flux density of optically thick blackbody
def blackbody2(fGHz, Tb, rAU=10 ) :
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
  
# --- gaussian function used to fit radial brightness --- #
def gaussian1d( x, a, b ) :
  return a * numpy.exp(-4.*math.log(2.) * numpy.power(x,2.)/pow(b,2))

# --- dump file of continuum and line emissivity (0-1) vs radius (in milliarcsec) at a particular frequency
# --- also calculate continuum and (peak) line flux, and Gaussian fit to FWHM of emission region
def radial (outfile, fGHz, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) :
    xAUarr = []
    eCarr = []
    eLarr = []
    sumC = 0.
    sumL = 0.
    ofile = open(outfile,"w")
    ofile.write("# alpha = %.2f\n" % alpha)
    ofile.write("# r0 = %.2f AU\n" % r0AU )
    ofile.write("# n0 = %.2e cm-3 at r0\n" % n0)
    ofile.write("# rmax = %.2e AU\n" % rmaxAU)
    ofile.write("# dx = %.4e AU\n" % dxAU)
    ofile.write("# f = %.2f GHz\n" % fGHz)
    ofile.write("# r (mas)   eC    eL\n")
    for xAU in numpy.arange(dxAU/2.,rmaxAU,dxAU) :
      EM = em(xAU, alpha, r0AU, n0, rmaxAU)
      tC = tauC(EM,fGHz)
      tL = tauL(EM,fGHz,fwhm)
      eC = 1.
      if (tC < 100.): 
        eC = (1. - math.exp(-1.*tC))
      sumC = sumC + eC * xAU * dxAU            # add (1-exp(-tauC)) * radius * dr
      eL = 1.
      if ((tL+tC) < 100.): 
        eL = 1. - math.exp(-1.*(tL+tC))
      sumL = sumL + eL * xAU * dxAU
      xAUarr.append( xAU )
      eCarr.append( eC )
      eLarr.append( eL )
    popt,pcov = curve_fit( gaussian1d, xAUarr, eCarr )
    for xAU,eC,eL in zip(xAUarr,eCarr,eLarr) :
      ofile.write("%8.2f  %8.5f  %8.5f  %8.5f\n" % (1000.*xAU/400., eC, gaussian1d(xAU, *popt),  eL))
    ofile.close()
    const = 4. * math.pi * kB * Te * 1.e44 / (dAU * dAU * c * c)
    return const*sumC*fGHz*fGHz, const*sumL*fGHz*fGHz, popt[1]
       # popt = amplitude and FWHM of gaussian fit to emissivity vs radius in AU

    # mas assumes distance of 400 pc

def radialModel( outfile, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) :
    fout = open( outfile, "w" ) 
    for freq in range(10,400,10) :
      fluxC,fluxL,FWHM_AU = radial("dummy", freq, fwhm, alpha, r0AU, n0, rmaxAU, dxAU )
      fout.write("%8.2f %8.3f  %8.2f %8.3f  %8.2f %8.3f  %8.2f\n" % (freq, math.log10(freq), fluxC, math.log10(fluxC), \
         fluxL, math.log10(fluxL), FWHM_AU*2.5) )    # 2.5 = 1000/400 converts AU into milliarcsec
    fout.close()
  
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
    for yAU in numpy.arange( dyAU/2., ymaxAU, dyAU ) :
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
    for logf in numpy.arange(0.4, 3.6, 0.1):
      fGHz = pow(10., logf)
      print fGHz
      logf = math.log10(fGHz)
      sumC = 0.
      sumL = 0.
      for xAU in numpy.arange(dxAU/2.,rmaxAU,dxAU) :
        EM = em(xAU, alpha, r0AU, n0, rmaxAU)    # emission measure in annulus at radius xAU

        tC = tauC(EM,fGHz)					   # continuum opacity for this emission measure
        eC = 1.
        if (tC < 100.): 
          eC = (1. - exp(-1.*tC))	               # 1-exp(-tauC) for this anner --- #
        sumC = sumC + eC * xAU * dxAU            # add (1-exp(-tauC)) * radius * dr

        tL = tauL(EM,fGHz,fwhm)		           # opacity at center of gaussian line
        eL = 1.
        if ((tL+tC) < 100.):
          eL = (1. - exp(-1.*(tL+tC)))
        sumL = sumL + eL * xAU * dxAU
      
      SC = const * sumC * fGHz*fGHz
      SL = const * sumL * fGHz*fGHz       
      ofile.write("%8.3f  %8.4f  %8.3f  %8.4f  %8.3f  %8.4f\n" %  \
         (fGHz,logf,SC,math.log10(SC),SL,math.log10(SL)))

# --- find gauss at velocity v relative to velocity 0, for width fwhm --- #
def gauss( v, fwhm ) :
    return exp(-2.7726 * pow(v/fwhm, 2.))

# freq for Halpha(n+1 -> n); for H and heavy element (K); velocity difference between them
def freq( n ):
    RH = 3.28805129e6	 # GHz
    RK = 3.28984202e6	 # GHz; for heavy atom like K
    fK = (RK * (1./float(n*n) - 1./float((n+1)*(n+1))))
    print fH, fK, -3.e5*(fK-fH)/fH

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
    S10AU = blackbody2( fGHz, Tb )			# flux density in Jy for optically thick emission 10 AU radius 
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


alphaList = numpy.arange(-2.,-4.,-.3)
r0List = numpy.arange(1.,10.,2.)
n0List = numpy.arange(1.e7,9.e7,2.e7)

def findBestModel( alphaList, r0List, n0List, fwhm_kms=30., rmaxAU=70., dxAU=0.5, fluxFile="../Figures/BNintFluxSummary", \
      sizeFile="../Figures/radialFit.dat" ) : 

  # read in measured fluxes and sizes
    b = numpy.loadtxt( fluxFile )     # freq, flux, unc
       # b[0] = freq, b[1] = intFlux, b[2] = fluxUncertainty
    c = numpy.loadtxt( sizeFile )
       # c[0] = freq, c[1] = FWHM, c[2] = FWHMuncertainty

  # now loop through all possible models, computing separate chisq for flux, size, and combined
    n = len(alphaList)*len(r0List)*len(n0List)
    print "evaluating %d models" % n
    chisqFWHM = numpy.zeros(n)
    chisqFlux = numpy.zeros(n)
    chisqBoth = numpy.zeros(n)
    j = 0
     
    for alpha in alphaList :
      for r0AU in r0List :
        for n0 in n0List : 
          j = j + 1
          for freq,flux,unc in zip(b[0],b[1],b[2]) :
            fluxC,fluxL,FWHM_AU = radial("dummy", freq, fwhm_kms, alpha, r0AU, n0, rmaxAU, dxAU )
            chisqFlux[j] = chisqFlux[j] + pow( (flux-fluxC)/unc, 2. )
          for freq,FWHM,unc in zip(c[0],c[1],c[2]) :
            fluxC,fluxL,FWHM_AU = radial("dummy", freq, fwhm_kms, alpha, r0AU, n0, rmaxAU, dxAU )
            chisqFWHM[j] = chisqFWHM[j] + pow( (FWHM-2.5*FWHM_AU)/unc, 2. )
          print "%5.2f  %4.1f  %.2e  %8.3f  %8.3f" % (alpha, r0AU, n0,
            chisqFlux[j]/math.sqrt(len(b[0])), chisqFWHM[j]/math.sqrt(len(c[0])))

# def examples () :
# radial (outfile, fGHz, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) 
# radial ("r43.dat", 43., 39., -10., 5., 1.e7, 50., 0.5 ) 
# flux ( outfile,     fwhm, alpha,  r0AU,   n0,  rmaxAU,  dxAU ) :
# flux ("fnormal.dat", 26.,  -10.,  10.0,  2.e7,   20.,   0.5 )
# flux ("BNmodel.dat",      23.,  -5.,    9.5,  3.e7,   30.,   0.5 )   # save
# flux ("BNmodel1.dat",      23.,  -3.3,    7,  4.5e7,   30.,   0.5 )   # good fit to most of the spectrum
# flux ("HII.dat",     23.,  -20.,   9.5,  3.e7,   30.,   0.5 )
# flux ("Imodel.dat",   23.,  -10,   6.8,  8.5e7,   60.,   0.5 )   # save - good fit
# flux ("IHII.dat",   23.,  -100,   7.6,  7.1e7,   60.,   0.5 )   # save - good fit
# flux ("BNmodel2.dat",      23.,  -10.,    10.5,  2.8e7,   30.,   0.5 )   # good fit to 43,89,229 GHz points
# flux ("Imodel.dat",   23.,  -10,   6.8,  8.5e7,   60.,   0.5 )   # save - good fit
# flux ("classicHII.dat",   23.,  -20,   6.8,  1.e7,   60.,   0.5 ) 
# flux ("classicWind.dat",   23.,  -2,   .1,  1.e10,   60.,   0.5 )
# flux ("classicHyper.dat",   23.,  -3,   3,  1.e7,   60.,   0.5 )

#  ---- feb 2011 ---
# flux ("Imodel.26feb2011.dat",   23.,  -20, 6.8,  3e8,  50., 0.5 )

#   ---- jul 2012 ----
# radial(outfile,           fGHz, fwhm, alpha, r0AU, n0, rmaxAU, dxAU ) 
# radial("Imodel.radial.dat", 662., 23., -50., 6.8,   3e8,  20,  0.1)
# flux ("Imodel.dat",   23.,  -20, 6.8,  3e8,  10., 0.1 )
# flux ("Imodelextreme.dat", 40., -50., 6.8, 6e8, 10., 0.1 )
# flux ("Imodel3e8.dat",   23.,  -50, 6.8,  3e8,  10., 0.1 )
# flux ("Imodel3e9.dat",   23.,  -50, 6.8,  3e9,  10., 0.1 )
# flux ("Imodel3e10.dat",   23.,  -50, 6.8,  3e10,  10., 0.1 )

# flux ( outfile,     fwhm, alpha,  r0AU,   n0,  rmaxAU,  dxAU ) :
# flux ("Imodel.dat",   30.,  -50,   7.5,  1.5e8,   20.,   0.1 )   # save - good fit
# flux( "BNmodel2.dat",  30.,  -3.5,  7.4,  5e7,   50.,   0.1 )	# good fit feb2011

#  ---- nov 2018 ----
# flux( "BNmodel2.dat",  30.,  -3.5,  7.4,  5.e7,   50.,   0.1 )	# good fit feb2011
# radial("BNmodel2.radial.43.dat", 43., 30., -3.5, 7.4, 5.e7,  50.,  0.1)
# radial("BNmodel2.radial.86.dat", 97., 30., -3.5, 7.4, 5.e7,  50.,  0.1)
# radial("BNmodel2.radial.230.dat", 230., 30., -3.5, 7.4, 5.e7,  50.,  0.1)
# radial("BNmodel2.radial.340.dat", 340., 30., -3.5, 7.4, 5.e7,  50.,  0.1)
# radialModel("BNmodel2.radialModel", 30., -3.5, 7.4, 5.e7, 50, 0.1 )

#   ---- dec 2018 ---
# findBestModel( alphaList, r0List, n0List ) 

# ============================ mar 2020 ================================

# -----------------------------------------------------------------------------------------------------------#
# findTe computes electron temperature using formula (1) in Balser et al. 2015 
# note: assumes optically thin emission, same beam filling factor for both line and continuum
# note: Balser et al. formula was valid only near 9 GHz, where correction factor a ~ 0.98; to make
#    this general, assume Te=8000., compute correction factor for actual freq, then solve
# -----------------------------------------------------------------------------------------------------------#
def findTe( fGHz, Tline, Tcont, deltaVkms, y=.08 ) :
    Te_assumed = 8000.
    a = calc_a( fGHz, Te_assumed )
    # print pow(7103. * pow(fGHz,1.1) * (Tcont/Tline) / (deltaVkms * (1+y)), 0.87)
        # this is the actual Balser formula, valid only near 9 GHz
    print pow(6985./a * pow(fGHz,1.1) * (Tcont/Tline) / (deltaVkms * (1+y)), 0.87)

    # validate using entries in Tables 1,2,3 of Balser et al for G034.133+0.471: 
    #   TL = 108.65, deltaVkms=25.29, y=(8.15*20.)/(108.65*25.29)=.0593, Tc=1064.6, Te=7655
    # findTe( 9., 108.65, 1064.6, 25.29, y=.0593)  gives Te=7624 (OK)

#findTe( 353.6, .890, .292, 32.64, y=0. )
#findTe( 231.9, .261, .220, 37.23, y=0. )
#findTe(  99.0, .040, .047, 34.67, y=0. )

# -----------------------------------------------------------------------------------------------------------#
# correction factor a, from Mezger 1967, A.2
# -----------------------------------------------------------------------------------------------------------#
def calc_a( fGHz, TK ):
    a = 0.366 * pow(fGHz, 0.1) * pow(TK, -0.15) * (math.log(4.995e-2/fGHz) + 1.5*math.log(TK))
    return a
  
# -----------------------------------------------------------------------------------------------------------#
# predict line/continuum ratio from Wilson eqn 14.29
# -----------------------------------------------------------------------------------------------------------#
def line_to_cont( fGHz, Te, deltaVkms ) :
    a = calc_a( fGHz, Te )
    print "a = %.3f" % a
    return (6.985e3/a) * pow(fGHz,1.1) * pow(Te,-1.15) / deltaVkms

# print line_to_cont( 85.69, 8000., 35. ) 
# print line_to_cont( 99.022, 8000., 35. ) 
# print line_to_cont( 231.9, 8000., 35. ) 
# print line_to_cont( 353.6, 8000., 35. ) 

# tauC3 is from Prozesky and Smits 2019, eqn 12 and 13
def tauC3( EM, fGHz ) :
    gamma = 0.5772
      # Euler-Mascheroni constant from Wikipedia
    gaunt = (math.sqrt(3)/math.pi) * math.log( pow(2.*kB*Te/(gamma*mass_e), 1.5) * mass_e/(math.pi*gamma*pow(e_esu,2)*fGHz*1.e9))
      # avg Gaunt factor from eqn 13
    tauC = (EM/pow(fGHz*1.e9,2)) * 3.08e18 * 8. * pow(e_esu,6) / (3 * math.sqrt(3) * pow(mass_e,3) * clight) * math.sqrt(math.pi/2) \
	* pow( mass_e/(kB*Te), 1.5) * gaunt
      # kappa from eqn 13, multiplied by length in cm (converted from pc)
    return gaunt,tauC

def tauMORELI( EM, fGHz ):
  # compute Gaunt factor using eqn A.1
    # gammasq = hplanck * (13.6/4.135e-6) * 1.e9 /(kB*Te)
      # this formula gives same result as one used below
    gammasq = 2.1798e-11 / (kB*Te)
      # Rydberg in ergs/kT
    print gammasq
    factor = hplanck*fGHz*1.e9/(kB*Te)
    x = 0.5 * factor * (1.+math.sqrt(10.*gammasq)) 
    a = 1.20 * math.exp(-1.*pow( (math.log(gammasq) - 1.)/3.7, 2) )
    b = 0.37 * math.exp(-1.*pow( (math.log(gammasq) - 1.)/2., 2) )
    gaunt = math.sqrt( pow( math.sqrt(3)/math.pi * math.exp(x) * kn(0,x), 2) \
      + pow( a - b*math.log(factor), 2) )
    #print gaunt, math.log10(gaunt)
    tauMORELI = 1.77e-12 * EM * 1.e12 * 3.08e16 * gaunt /(pow(fGHz*1.e9,2)*pow(Te,1.5))
       # formula given was in SI units
    return gaunt, tauMORELI
    
# Brocklehurst 1972
def tauBr( EM, fGHz ) :
    gamma = 0.5772
    factor1 = kB*Te/(mass_e * pow(clight,2))
    re = pow(e_esu,2)/(mass_e * pow(clight,2))
    gaunt = math.log( pow(2/gamma, 2.5) * pow(factor1, 1.5) * clight/(fGHz * 1.e9  * 2. * math.pi * re))
    tauBr = EM * pc * pow(re,3) * pow(clight/(fGHz*1.e9),2) * (8./(3.*math.sqrt(2.*math.pi))) * pow(factor1,-1.5) * gaunt
    return gaunt, tauBr
    
# compare eqns 10.33/10.34 of Rolhfs with eqn 10.35
def cmpTauC( EM, fGHz) :

  gff = math.log(4.955e-2/fGHz) + 1.5*math.log(Te)
  tauC1 = 3.014e-2 * pow(Te,-1.5) * pow(fGHz, -2) * EM * gff
  print "\nRohlfs eqn 10.33,10.34:"
  print "   gaunt = %.4f" % gff
  print "   tauC = %.3e" % tauC1

  print "\nBrocklehurst 1972:"
  gaunt, tau6 = tauBr( EM, fGHz )
  print "   gaunt = %.4f" % gff
  print "   tauC = %.3e" % tau6

  print "\nRohlfs eqn 10.35:"
  print "   tauC = %.3e" % ( tauC(EM,fGHz) )

  print "\nProzesky and Smits 2019, eqn. 12,13:"
  gaunt, tauC2 = tauC3( EM, fGHz )
  print "   gaunt = %.4f" % gaunt
  print "   tauC = %.3e" % tauC2

  print "\nMORELI, eqn A.1, 3"
  gaunt, tauCM = tauMORELI( EM, fGHz)
  print "   gaunt = %.4f" % gaunt
  print "   tauC = %.3e" % tauCM
  
# cmpTauC( 1.e4, 200. )

# pressure broadening from table B.1 of Baez-Rubio 2013
def delnu1(n, Ne, Te) :
    delnu1_BS = 1.9e-8 * pow(n,4.4) * Ne * pow(Te,-0.1)  # n > 100
    delnu1_W = 6.7e-9 * pow(n,4.6) * Ne * pow(Te,-0.1)
    delnu1_S = 8.e-10 * pow(n,5) * Ne * pow(Te,-0.1)
    return delnu1_BS,delnu1_W,delnu1_S

'''
for Ne in [ 1.e4, 1.e5, 1.e6, 1.e7, 1.e8] :
  print " "
  for nnn,fGHz in zip( [26,30,40,52],[353.,232.,99.,88.]) :
    d1,d2,d3 = delnu1( nnn, Ne, 8000.) 
    dv1 = d1/(fGHz*1.e9) * 3.e5
    dv2 = d2/(fGHz*1.e9) * 3.e5
    dv3 = d3/(fGHz*1.e9) * 3.e5
    print "%.2e  %d  %6.2f  %6.2f  %6.2f" % (Ne,nnn,dv1,dv2,dv3)
'''
# ===== oct 2020 =====

# -----------------------------------------------------------------------------------------------------------#
# model dictionary provides the parameters to compute the emission measure projected on the (x,y) plane
#       for a biconical source with axis of neutral disk exactly normal to the line of sight; also includes
#       electron temperature and fwhm of recomb lines for convenience (not needed to compute em),
#       and distance of source in AU
#
#   alpha = radial dependence of electron density, r^alpha
#   r_core_AU = radius of core; electron density is constant inside this, drops as r^alpha outside
#   ne_core =  electron density for r <= r_core_AU
#   phi_deg = half-opening polar angle of ionized cones; phi_deg = 90 for spherical outflow
#   incr_AU = the step size in x,y,z
#   rmax_AU = the max radius in AU of the calculation (note that incremental elements are included only
#      if their radius is less than rmax_AU, so entire x,y,z range is not included; if rmax_AU is too small,
#      low freq flux densities are underestimated and flux curve appears to have an inflection point;
#      for final fits, be sure to use rmax_AU of 100 AU
#   Te = electron temperature in K (included in model for convenience only; not needed to compute EMimg)
#   fwhm = FWHM in km/sec of Halpha recombination lines (included in model for convenience only)
#   dist_AU = distance of source in AU (assume 400 AU for Orion)
# -----------------------------------------------------------------------------------------------------------#

# this is what I published in 2013
modelA = { "alpha" :  -3.5,
           "r_core_AU" : 7.4,
           "ne_core" : 5.e7,
           "phi_deg" : 90.,
	   "incr_AU" : 1.,
           "rmax_AU" : 100.,
           "Te" : 8000.,
           "fwhm" : 35.,
           "dist_AU" : 400.,
         }  

modelB = { "alpha" :  -3.5,
           "r_core_AU" : 8.2,
           "ne_core" : 9.e7,
           "phi_deg" : 90.,
	   "incr_AU" : 1.,
           "rmax_AU" : 50.,
           "Te" : 8000.,
           "fwhm" : 35.,
           "dist_AU" : 400.,
         }  

modelC = { "alpha" :  -2.6,
           "r_core_AU" : 3.8,
           "ne_core" : 2.e8,
           "phi_deg" : 75.,
	   "incr_AU" : 1.,
           "rmax_AU" : 50.,
           "Te" : 8000.,
           "fwhm" : 35.,
           "dist_AU" : 400.,
         }  

# modelD gives spectral slope of 1.2 that is observed
# to get greater spread between recomb peak and continuum, increase ne
modelD = { "alpha" :  -3.125,
           "r_core_AU" : 5.,
           "ne_core" : 1.e8,
           "phi_deg" : 90.,
	   "incr_AU" : 1.,
           "rmax_AU" : 100.,
           "Te" : 8000.,
           "fwhm" : 35.,
           "dist_AU" : 400.,
         }  

# modelE is pure wright and barlow fit
modelE = { "alpha" :  -3.125,
           "r_core_AU" : 2.3,
           "ne_core" : 1.e9,
           "phi_deg" : 90.,
	   "incr_AU" : 1.,
           "rmax_AU" : 100.,
           "Te" : 8000.,
           "fwhm" : 35.,
           "dist_AU" : 400.,
         }  

# modelF is a lower density, more extended source to fit recomb line peaks
modelF = { "alpha" :  -2.8,
           "r_core_AU" : 6.,
           "ne_core" : 3.e7,
           "phi_deg" : 90.,
	   "incr_AU" : 1.,
           "rmax_AU" : 60.,
           "Te" : 8000.,
           "fwhm" : 35.,
           "dist_AU" : 400.,
         }  

# modelG is optically thin spherical HII region
modelG = { "alpha" :  -10.,
           "r_core_AU" : 10.,
           "ne_core" : 1.e6,
           "phi_deg" : 90.,
	   "incr_AU" : 1.,
           "rmax_AU" : 20.,
           "Te" : 8000.,
           "fwhm" : 35.,
           "dist_AU" : 400.,
         }  

# -----------------------------------------------------------------------------------------------------------#
# EMcalc computes projected emission measure toward one quadrant of the source (x>0, y>0)
# creates new pickle file, writes model and EMing to it
# -----------------------------------------------------------------------------------------------------------#
def EMcalc( model, pklfile ):
    phi_rad = numpy.radians( model["phi_deg"] )
    xAU = yAU = zAU = numpy.arange( model["incr_AU"]/2., model["rmax_AU"] + model["incr_AU"], model["incr_AU"] )
    nn = len( xAU )
    EMimg = numpy.zeros( (nn,nn) )
    for i in range( 0, nn ):
      for j in range( 0, nn ):
        for k in range( 0, nn ):
          r = math.sqrt( pow(xAU[i],2) + pow(yAU[j],2) + pow(zAU[k],2) ) 
          if r > model["rmax_AU"] :
            break
        # changed to make sure core is filled between the bicones
          # elif (r > model["r_core_AU"]) and (math.atan( math.sqrt(pow(xAU[i],2) + pow(zAU[k],2))/yAU[j] ) > phi_rad) :
          elif (math.atan( math.sqrt(pow(xAU[i],2) + pow(zAU[k],2))/yAU[j] ) > phi_rad) :
            break
          else :
            if r < model[ "r_core_AU" ] :
              ne = model[ "ne_core" ]
            else:
              ne = model[ "ne_core" ]  * pow(r/model[ "r_core_AU"], model["alpha"] )
            EMimg[i,j] = EMimg[i,j] + 2.* ne * ne * model["incr_AU"] * 1.496e13 / 3.0856e18
		# add emission measure increment for this cell, in cm-6 * pc
                # factor of 2 accounts for the contribution from the mirror image cell at z < 0

  # append EM image file to model dictionary and pickle it
    model["EMimg"] =  EMimg   
    fout = open( pklfile, "wb" )
    pickle.dump( model, fout )
    fout.close()

# -----------------------------------------------------------------------------------------------------------#
# mirror image from one quadrant to 4 quadrants (for now, requires image to be square)
# -----------------------------------------------------------------------------------------------------------#
def expandImage( img1 ) :
    nn = len(img1[0])
    img2 = numpy.zeros( (2*nn, 2*nn) )
    for i in range(0, nn ):
      for j in range(0, nn ):
         img2[nn+i,nn+j] = img1[i,j]       # e.g, for nn=10, img2[10,10] = img1[0,0]
         img2[nn-i-1,nn+j] = img1[i,j]     # e.g, img2[8,10] = img1[1,0]
         img2[nn-i-1,nn-j-1] = img1[i,j]   
         img2[nn+i,nn-j-1] = img1[i,j]
    return img2
   
# -----------------------------------------------------------------------------------------------------------#
# given EM image (usually just 1 quadrant), return continuum and Ha recomb line brightness temperature images
#  fGHz = single frequency in GHz
#  Te = electron temperature in K
#  fwhm = linewidth of recombination lines in km/sec
#  EMimg = emission measure image (units cm^-6 * pc) for all x,y
# -----------------------------------------------------------------------------------------------------------#
def Tff ( fGHz, Te, fwhm, EMimg ) :
    dfkHz = 1.e11 * fwhm * fGHz / clight
    factor1 =  3.014e-2 * pow(Te, -1.5) * pow(fGHz, -2.) * gaunt(fGHz, Te) 
    factor2 = 1.92e3 * pow(Te,-2.5) / dfkHz
    TcontImage = Te * ( 1. - numpy.exp(-1.* factor1 * EMimg) )
    TlineImage = Te * ( 1. - numpy.exp(-1.* (factor1 + factor2) * EMimg) )
    return TcontImage, TlineImage

# -----------------------------------------------------------------------------------------------------------#
# 2D gaussian fit to synthetic image
# copied from https://scipython.com/blog/non-linear-least-squares-fitting-of-a-two-dimensional-data/
# -----------------------------------------------------------------------------------------------------------#
def gauss2d( x, y, zmax, fwhm_x, fwhm_y ) :
    a = -4.*math.log(2.)/pow(fwhm_x,2)
    b = -4.*math.log(2.)/pow(fwhm_y,2)
    return zmax * numpy.exp( a*numpy.power(x,2) + b*numpy.power(y,2) )

def _gaussian( M, *args ):
    x,y = M
    arr = gauss2d(x,y,*args[0:3])
    return arr

def fit2dgaussian(x,y,img):
    xdata = numpy.vstack( (x.ravel(), y.ravel()) )    # turn into 1-d problem
    ncent = int(len(img[0])/2 + 0.5)
    guesspeak = img[ncent,ncent]
    popt, pcov = curve_fit( _gaussian, xdata, img.ravel(), [guesspeak,.06,.04] )
    #print "peak = %.1f mJy" % popt[0]
    #print "maj = %.3f arcsec" % popt[1]
    #print "min = %.3f arcsec" % popt[2]
    return popt

# -----------------------------------------------------------------------------------------------------------#
# pre-convolve with gaussian, fit 2d gaussian to result, return maj,min corrected for smoothing gaussian
# this routine is needed because directly fitting gaussian to bicone gives nonsensical answwers
# -----------------------------------------------------------------------------------------------------------#
def fitConv2dgaussian(x, y, img1):
    conv = .02
    gaussBeam = gauss2d( x, y, 1., conv, conv)
    img2 = convolve2d( img1, gaussBeam, mode='same' )
    pfit = fit2dgaussian( x, y, img2 )
    bmaj = math.sqrt( pow(pfit[1],2) - pow(conv,2) )
    bmin = math.sqrt( pow(pfit[2],2) - pow(conv,2) )
    #print "fit bmaj = %.3f --> %.3f" % (pfit[1],bmaj)
    #print "fit bmin = %.3f --> %.3f" % (pfit[2],bmin)
    print "bmaj, bmin = %.3f, %.3f" % (bmaj,bmin)
    return bmaj,bmin,img2

# -----------------------------------------------------------------------------------------------------------#
# test gaussian fitting algorithm
# fit2dgaussian fails when fitting unsmoothed biconical source, so I convolve the model image with a 
#   round gaussian, fit a gaussian to the resultant, then correct sizes for convolving beam
# this routine is designed to test this procedure
# -----------------------------------------------------------------------------------------------------------#
def testgaussfit( bmaj, bmin ) :
    dx = bmin/5.
    xmax = 5.*bmaj
    xp = numpy.arange( dx/2., xmax + dx/2., dx )   # positive x values
    x_arcsec = numpy.concatenate( (-1.*xp[::-1],xp) )  
    x,y = numpy.meshgrid( x_arcsec, x_arcsec )
    img1 = gauss2d( x, y, 1., bmaj, bmin)
    bmajfit,bminfit = fitConv2dgaussian( x, y, img1 )
    print "bmaj = %.3f, bmajfit = %.3f" % (bmaj,bmajfit)
    print "bmin = %.3f, bminfit = %.3f" % (bmin,bminfit)

# -----------------------------------------------------------------------------------------------------------#
# given a model containing an emission measure file, plot brightness temperature image 
#   at freq fGHz; if type='continuum', plot continuum brightness; if type='line', plot Halpha;
#   in a second panel, plot actual image
# -----------------------------------------------------------------------------------------------------------#
def plotTb (pklfile, fGHz, type='continuum', rotdeg=70. ) :
    try :
      fin = open( pklfile, "rb" )
      model = pickle.load( fin )
      fin.close()
    except :
      print "FATAL: couldn't find pickle file %s" % pklfile

  # compute continuum and line brightness temp images for frequency fGHz
    TcontImg, TlineImg = Tff ( fGHz, model["Te"], model["fwhm"], model["EMimg"] ) 

  # mirror to 4 quadrants, create x and y axes in arcsec
    if type == 'line' :
      img1 = expandImage (TlineImg)
    else :
      img1 = expandImage (TcontImg)
    xAU = yAU = numpy.arange( model["incr_AU"]/2., model["rmax_AU"] + model["incr_AU"], model["incr_AU"] )
    x_arcsec = numpy.concatenate( (-1. * xAU[::-1],xAU) )/model["dist_AU"]
    y_arcsec = numpy.concatenate( (-1. * yAU[::-1],yAU) )/model["dist_AU"]
    x,y = numpy.meshgrid( x_arcsec, y_arcsec )
  
  # gaussian smooth by 2 cells (not found to be useful)
  # cell = 2.*model["incr_AU"]/model["dist_AU"]
  # print "convolve raw image with %.3f arcsec gaussian to smooth" % cell
  # imgsmooth = gauss2d( x, y, 1., cell, cell ) 
  # img2 = convolve2d( img1, imgsmooth, mode='same' )

  # fitting Gaussian size of bicone (first convolve with gaussian, then fit, then deconvolve - fails otherwise)
    bmaj,bmin,img2 = fitConv2dgaussian(x, y, img1)
    img3 = gauss2d( x, y, 1., bmaj, bmin ) 

  # rotate both model image and gaussian fit by rotdeg
    img4 = rotate(img1, rotdeg, reshape=False )    # rotate model image
    img5 = rotate(img3, rotdeg, reshape=False )    # rotate Gaussian fit
    # print "size of img1, img2 = ", img1.shape, img2.shape
 
    fig = pyplot.figure()
    ax1 = fig.add_axes( [.1,.1,.7,.7] )
    arcsecMax = model["rmax_AU"]/model["dist_AU"]   # convert from AU to arcsec
    #ax1.axis ( [-arcsecMax,arcsecMax,-arcsecMax,arcsecMax] )
    print "checking plot range: ",-arcsecMax,arcsecMax,x_arcsec[0],x_arcsec[-1]
    ax1.axis( [x_arcsec[0],x_arcsec[-1],y_arcsec[0],y_arcsec[-1]] )
    cs = ax1.imshow( img4, origin='lower', extent=[-arcsecMax,arcsecMax,-arcsecMax,arcsecMax] )
    #ax1.contour( x, y, img5, colors='white'  )
    ax1.contour( x, y, img5, [0.5], colors='white'  )
    fig.colorbar( cs )
    ax1.annotate( "%.0f GHz" % fGHz, (.04,.95), xycoords='axes fraction', color='white')
    ax1.annotate( "%s" % type, (.04,.91), xycoords='axes fraction', color='white')
    
    pyplot.show()
       
# -----------------------------------------------------------------------------------------------------------#
# given a model containing an emission measure file, compute flux at frequencies in fGHzArray
# -----------------------------------------------------------------------------------------------------------#
def computeFlux( pklfile, fGHzArray, outFile=None ):
    try :
      fin = open( pklfile, "rb" )
      model = pickle.load( fin )
      fin.close()
    except :
      print "FATAL: couldn't find pickle file %s" % pklfile

  # pixel size in arcsec is incr_AU/400. x incr_AU/400.
    dOmega_per_pixel = pow( model["incr_AU"]/dist_AU * math.pi/(60.*60.*180.), 2 )
  # print "dOmega_per_pixel = %.3e" % dOmega_per_pixel

    ContinuumFluxArray = []
    HalphaPeakFluxArray = []
    for fGHz in fGHzArray:

    # convert from emission measure to brightness temperature
      TcontImg, TlineImg = Tff ( fGHz, model["Te"], model["fwhm"], model["EMimg"] ) 

    # total flux = 2kT fHz^2/c^2 * dOmega * 10^-23 Jy/(ergs s^-1 cm^-2 Hz^-1) summed over image 
      factor = 2. * kB * pow( (1.e9*fGHz)/clight, 2) * dOmega_per_pixel / 1.e-23
      #print "factor = %.3e" % factor
      continuum_flux_mJy = 4. * 1.e3 * factor * numpy.sum( TcontImg )
      Halpha_peak_line_flux_mJy = 4. * 1.e3 * factor * numpy.sum( TlineImg )
	# multiply by 4 because image covers just 1 quadrant of source
	# multiply by 1000 to convert to mJy
      # print fGHz, continuum_flux_mJy, Halpha_peak_line_flux_mJy
      ContinuumFluxArray.append( continuum_flux_mJy )
      HalphaPeakFluxArray.append( Halpha_peak_line_flux_mJy )

  # rewrite the pickle file, including new dictionary items 
    model[ "fGHzArray" ] = fGHzArray
    model[ "contFlux_mJy" ] = numpy.array( ContinuumFluxArray )
    model[ "lineFlux_mJy" ] = numpy.array( HalphaPeakFluxArray )
    fout = open( pklfile, "wb" )
    pickle.dump( model, fout )
    fout.close()

  # optional: write fluxes and line/cont ratios to text file if filename is given
    if (outFile) :
      fout = open( outFile, "w" )
      fout.write("# %s\n" % outFile)
      fout.write("#   alpha = %.2f\n" % model["alpha"] )
      fout.write("#   r_core_AU = %.2f\n" % model["r_core_AU"] ) 
      fout.write("#   ne_core = %.2e\n" % model["ne_core"] )
      fout.write("#   phi_deg = %.2e\n" % model["phi_deg"] )
      fout.write("#   Te = %.0f\n" % model["Te"] )
      fout.write("#   fwhm = %.1f\n" % model["fwhm"] )
      fout.write("#   dist_AU = %.0f\n" % model["dist_AU"] )
      fout.write("#     fGHz     ContFlux_mJy  AlphaFlux_mJy  line/cont\n#\n")
      for fGHz, ContFlux, HalphaFlux in zip(fGHzArray,ContinuumFluxArray,HalphaPeakFluxArray) :
        fout.write("%10.3f  %10.3f  %12.3f  %12.3f\n" % (fGHz,ContFlux,HalphaFlux,(HalphaFlux-ContFlux)/ContFlux))
      fout.close()
    fout.close()
       
# -----------------------------------------------------------------------------------------------------------#
# overplot data on theoretical flux curves (for continuum and peak Halpha)
# also compute chisq, deviation of measured from theoretical
# -----------------------------------------------------------------------------------------------------------#

def plotFluxes( pklfile, measFluxFile ) :
  # read model fluxes
    try :
      fin = open( pklfile, "rb" )
      model = pickle.load( fin )
      fin.close()
    except :
      print "FATAL: couldn't find pickle file %s" % pklfile

  # read data from csv file (downloaded from Google Drive)
    data = numpy.genfromtxt( measFluxFile, dtype=None, delimiter=",", names=True)
    contdata = numpy.genfromtxt( "BNcontFlux.csv", dtype=None, delimiter=",", names=True)

  # open the figure
    pp = PdfPages("BNfluxes.pdf")
    fig = pyplot.figure( figsize=(11,8) )
    #fig = pyplot.figure( figsize=(8,8) )
    ax = fig.add_subplot( 111 )
    ax.set_xlim( 3., 800. )
    ax.set_ylim( 1., 2000. )

  # plot models first so points will lie above them
    ax.plot( model["fGHzArray"], model["contFlux_mJy"] )
    ax.plot( model["fGHzArray"], model["lineFlux_mJy"], '--' )

  # fit polynomial (linear works best) to the measured continuum data
    #pcoeff = numpy.polyfit( numpy.log10(contdata["f_GHz"]), numpy.log10(contdata["flux_mJy"]), deg=1 )
    #print "polynomial fit to observed continuum fluxes : ",pcoeff
    #ax.plot( model["fGHzArray"], numpy.power(numpy.polyval(pcoeff, numpy.log10(model["fGHzArray"])),10.), color='red' )
        # why does this fail??

    pcoeff = numpy.polyfit( numpy.log(contdata["f_GHz"]), numpy.log(contdata["flux_mJy"]), deg=1 )
    print "polynomial fit to observed continuum fluxes : ",pcoeff
    #ax.plot( model["fGHzArray"], numpy.exp(numpy.polyval(pcoeff, numpy.log(model["fGHzArray"]))), color='red' )
        # but this works fine!
    
  # create cubic spline fit to the model -- use to get fluxes at exactly the measured freqs
    fcont = interp1d( model["fGHzArray"], model["contFlux_mJy"] )
    fline = interp1d( model["fGHzArray"], model["lineFlux_mJy"] )
    chisqCont = chisqAlpha = 0.
    nCont = nAlpha = 0

  # plot points one at a time because errorbar doesn't allow array of symbols or colors
    for n in range(0, len(data)) :
      ax.errorbar( data["f_GHz"][n], data["flux_mJy"][n], yerr=data["unc"][n], 
        marker=data["marker"][n], markersize=data["markersize"][n], color=data["color"][n],
        markeredgecolor='k', elinewidth=1.5, capthick=1.5, alpha=0.7 )

    # annotate recomb lines, indicated by alpha or beta in data["ref"]
      if "alpha" in data["text"][n] or "beta" in data["text"][n]  :
        # print data["text"][n]
        ax.annotate( data["text"][n], (data["f_GHz"][n],data["flux_mJy"][n]), \
           xytext=(data["xoff"][n],data["yoff"][n]), textcoords='offset points', 
           color=data["color"][n] )

    # add to chisqCont or chisqAlpha (ignore beta lines)
      frq = data["f_GHz"][n]
      flx = data["flux_mJy"][n]
      if "alpha" in data["text"][n] :
        #print "   %7.3f  %7.2f  %7.2f  %7.2f  %6.2f" % (frq, flx, fline(frq), flx-fline(frq), data["unc"][n])
        chisqAlpha = chisqAlpha + pow( (flx - fline(frq)), 2 )    # data["unc"] = 0 for lines
        nAlpha = nAlpha + 1
      elif "beta" in data["text"][n] :
        continue
      else :
        #print "   %7.3f  %7.2f  %7.2f  %7.2f  %6.2f" % (frq, flx, fcont(frq), flx-fcont(frq), data["unc"][n])
        chisqCont = chisqCont + pow( (flx - fcont(frq))/data["unc"][n], 2 )
        nCont = nCont + 1
        
    # this would be the nice way to compute chisq if all points were continuum
    # chisqCont = numpy.sum( numpy.power( (data["f_GHz"]-fcont(data["f_GHz"]))/data["unc"], 2 ) )/ len(data["f_GHz"])
    # ax.plot( data["f_GHz"], fcont(data["f_GHz"]), "o" )

    print "chisqCont  = %.0f" % (chisqCont/nCont)
    print "chisqAlpha = %.0f" % (chisqAlpha/nAlpha)

    ax.annotate( "alpha = %.2f" % model["alpha"], (.04,.95), xycoords='axes fraction')
    ax.annotate( "r_core_AU = %.2f" % model["r_core_AU"], (.04,.91), xycoords='axes fraction')
    ax.annotate( "ne_core = %.1e" % model["ne_core"], (.04,.87), xycoords='axes fraction')
    ax.annotate( "phi_deg = %.1f" % model["phi_deg"], (.04,.83), xycoords='axes fraction')
    ax.annotate( "chisqCont = %.0f" % (chisqCont/nCont), (.04,.77), xycoords='axes fraction')
    ax.annotate( "chisqAlpha = %.0f" % (chisqAlpha/nAlpha), (.04,.73), xycoords='axes fraction')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel( 'frequency (GHz)', size=14 )
    ax.set_ylabel( 'integrated flux density (mJy)', size=14 )
    ax.tick_params( which='both', length=6, labelsize=14 )
    pyplot.grid(True)

    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()  


# -----------------------------------------------------------------------------------------------------------#
# plot measured/avg flux vs epoch to search for evidence of time variability
# -----------------------------------------------------------------------------------------------------------#
def timePlot( ) :
    contdata = numpy.genfromtxt( "BNcontFlux.csv", dtype=None, delimiter=",", names=True)

  # fit polynomial (linear works best) to the measured continuum data
    pcoeff = numpy.polyfit( numpy.log(contdata["f_GHz"]), numpy.log(contdata["flux_mJy"]), deg=1 )
    print "polynomial fit to observed continuum fluxes : ",pcoeff
  
  # compute fractional deviation from average plot
    meanFlux = numpy.exp(numpy.polyval(pcoeff, numpy.log(contdata["f_GHz"])))
    fracDev = (contdata["flux_mJy"] - meanFlux)/meanFlux
    fracError = contdata["unc"]/meanFlux

  # open the figure
    pp = PdfPages("timePlot.pdf")
    fig = pyplot.figure( figsize=(11,8) )
    ax = fig.add_subplot( 111 )
    ax.set_xlim( 1990, 2022 )
    ax.set_ylim( -.6,.6 )
    for i in range(0,len(contdata)) :
      if contdata["f_GHz"][i] > 30. :
        if contdata["f_GHz"][i] > 80. : 
          contdata["color"][i] = "r" 
        elif (contdata["f_GHz"][i] > 30.) and (contdata["f_GHz"][i] < 80.) : 
          contdata["color"][i] = "b"
        print contdata["f_GHz"][i], meanFlux[i], contdata["flux_mJy"][i], fracDev[i], fracError[i], \
          contdata["epoch"][i], contdata["color"][i]
        ax.errorbar( contdata["epoch"][i], fracDev[i], yerr=fracError[i], 
          marker=contdata["marker"][i], markersize=contdata["markersize"][i], color=contdata["color"][i],
          markeredgecolor='k', elinewidth=1.5, capthick=1.5, alpha=0.7 )


    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()  
    

    
  

fGHzArray = numpy.power(10., numpy.arange(0.2,3.8,0.1) )

#EMcalc( modelG, "modelG.pkl" )
#plotTb( "modelG.pkl", 43., type='continuum', rotdeg=0. )
#computeFlux( "modelG.pkl", fGHzArray, outFile="Flux_modelG.dat"  )
#plotFluxes( "modelG.pkl", "../FluxPlot/BN - IntegratedFluxes.csv" )

#plotTb( "modelE.pkl", 340., type='continuum', rotdeg=0. )
#plotTb( "modelE.pkl", 340., type='line', rotdeg=0. )
#plotTb( "modelC.pkl", 99., 0. )
#plotTb( "modelC.pkl", 220., 0. )
#plotTb( "modelC.pkl", 340., 0. )
#plotFluxes( "modelE.pkl", "../FluxPlot/BN - IntegratedFluxes.csv" )
#testgaussfit( .15, .04 )

#timePlot()

#fGHzArray = [85.69,99.022,231.9,353.6]
# computeFlux( "modelE.pkl", fGHzArray, outFile="Flux_modelE.dat"  )
