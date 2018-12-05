# BNmodel.py
# compute flux density and size for HII region near BN

import numpy
import math
from scipy.optimize import curve_fit
import cPickle as pickle

# ------------- constants ------------------ #
Te = 8000            # electron temp in K
kB = 1.3806e-16       # boltzmann's constant, erg K-1
hplanck = 6.626e-27        # plank's constant, erg s
clight = 2.9979e10   # c in cm/sec
au = 1.496e13        # 1 AU in cm
pc = 3.0856e18		 # 1 pc in cm
G = 6.6732e-8        # grav constant in dynes cm^2 g^-2
dAU = 400 * pc/au    # 400 pc in AU 
mass_e = 9.109e-28	 # electron mass in grams
mass_H = 1.6605e-24  # atomic mass unit in grams
mass_sun = 1.989e33  # solar mass in g

# --- gaussian function used to fit radial brightness --- #
def gaussian( x, a, b ) :
  return a * numpy.exp(-4.*math.log(2.) * numpy.power(x,2.)/pow(b,2))

# --- write out file of continuum and line emissivity (0-1) vs radius (in milliarcsec) at a particular frequency
# --- also calculate continuum and (peak) line flux, and Gaussian fit to FWHM of emission region
def radial (outfile, fGHz, fwhm_kms, alpha, r0AU, n0, rmaxAU, dxAU ) :
    xAUarr = []
    eCarr = []
    eLarr = []
    sumC = 0.
    sumL = 0.

  # compute radial profiles of emissivities
    for xAU in numpy.arange(dxAU/2.,rmaxAU, dxAU) :
      EM = em(xAU, alpha, r0AU, n0, rmaxAU)
      tC = tauC(EM,fGHz)
      tL = tauL(EM,fGHz,fwhm_kms)
      eC = 1.
      if (tC < 100.): 
        eC = (1. - math.exp(-1.*tC))
      sumC = sumC + eC * 2. * math.pi * xAU * dxAU            # add (1-exp(-tauC)) * 2 * pi * radius * dr
      eL = 1.
      if ((tL+tC) < 100.): 
        eL = 1. - math.exp(-1.*(tL+tC))
      sumL = sumL + eL * 2. * math.pi * xAU * dxAU
      xAUarr.append( xAU )
      eCarr.append( eC )
      eLarr.append( eL )

  # do gaussian fit to emissivity profile to get FWHM
    popt,pcov = curve_fit( gaussian, xAUarr, eCarr )
       # popt = amplitude and FWHM of gaussian fit to emissivity vs radius in AU
    FWHM_milliarcsec = 2.5*popt[1]
       # 400 AU = 1000 milliarcsec

  # write out emissivity (0-1) vs radius in AU
    ofile = open(outfile,"w")
    ofile.write("# alpha = %.2f\n" % alpha)
    ofile.write("# r0 = %.2f AU\n" % r0AU )
    ofile.write("# n0 = %.2e cm-3 at r0\n" % n0)
    ofile.write("# rmax = %.2e AU\n" % rmaxAU)
    ofile.write("# dx = %.4e AU\n" % dxAU)
    ofile.write("# f = %.2f GHz\n" % fGHz)
    ofile.write("# r (mas)   eC    eL\n")
    for xAU,eC,eL in zip(xAUarr,eCarr,eLarr) :
      ofile.write("%8.2f  %8.5f  %8.5f  %8.5f\n" % (1000.*xAU/400., eC, gaussian(xAU, *popt),  eL))
    ofile.close()

  # compute integrated flux density in mJy = 2kT/lambda^2 * d(Omega)
  # d(Omega) = (area of src in AU^2)/(distance to src in AU)^2 in steradians
  # lambda = clight/(1.e9*fGHz) in cm 
  # convert from cgs units to Jy: x 1.e23
  # convert from Jy to mJy: x 1.e3
  # old conversion factor: const = 2 * kB * Te * 1.e44 / (dAU * dAU * c * c)
    S_continuum = 2 * kB * Te * sumC/(dAU*dAU) * pow(fGHz*1.e9/clight,2) * 1.e26
    S_linePlusContinuum = 2 * kB * Te * sumL/(dAU*dAU) * pow(fGHz*1.e9/clight,2) * 1.e26
    return S_continuum, S_linePlusContinuum, FWHM_milliarcsec

def radialModel( outfile, fwhm_kms, alpha, r0AU, n0, rmaxAU, dxAU ) :
  fout = open( outfile, "w" ) 
  for freq in range(10,400,10) :
    fluxC,fluxL,FWHM = radial("dummy", freq, fwhm_kms, alpha, r0AU, n0, rmaxAU, dxAU )
    fout.write("%8.2f %8.3f  %8.2f %8.3f  %8.2f %8.3f  %8.2f\n" % (freq, math.log10(freq), fluxC, math.log10(fluxC), \
       fluxL, math.log10(fluxL), FWHM) )   
  fout.close()
  
# --- continuum opacity, Rohlfs & Wilson 9.35 --- #
def tauC( EM, fGHz ) :
  return (8.235e-2 * pow(Te, -1.35) * pow(fGHz, -2.1) * EM)

# --- line opacity, Rohlfs & Wilson 13.27 --- #
def tauL( EM, fGHz, fwhm ) :
  dfkHz = 1.e11 * fwhm * fGHz / clight
  return (1.92e3 * pow(Te, -2.5) * EM / dfkHz )

# --- find emission measure at projected distance xAU from the center --- #
def em (xAU, alpha, r0AU, n0, rmaxAU) :
  dyAU = .1
  ymaxAU = math.sqrt(rmaxAU*rmaxAU - xAU*xAU)		# integrate to same max radius each time
  em = 0.
  for yAU in numpy.arange( dyAU/2., ymaxAU, dyAU ) :
    rAU = math.sqrt(xAU*xAU + yAU*yAU)
    ne = n0
    if (rAU > r0AU) :
      ne = n0 * pow(rAU/r0AU, alpha)	# electron density at this radius, cm-3
    em = em + ne * ne * dyAU * au/pc    # emission measure in cm-6 pc
  return 2.*em							# x 2 because we integrated through half the column  

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
    params = numpy.zeros( (n,3) )
    j = 0
     
    for alpha in alphaList :
      for r0AU in r0List :
        for n0 in n0List : 
          params[j] = [alpha,r0AU,n0]
          for freq,flux,unc in zip(b[0],b[1],b[2]) :
            fluxC,fluxL,FWHM_AU = radial("dummy", freq, fwhm_kms, alpha, r0AU, n0, rmaxAU, dxAU )
            chisqFlux[j] = chisqFlux[j] + pow( (flux-fluxC)/unc, 2. )
          chisqFlux[j] = chisqFlux[j]/len(b[0])
          for freq,FWHM,unc in zip(c[0],c[1],c[2]) :
            fluxC,fluxL,FWHM_AU = radial("dummy", freq, fwhm_kms, alpha, r0AU, n0, rmaxAU, dxAU )
            chisqFWHM[j] = chisqFWHM[j] + pow( (FWHM-2.5*FWHM_AU)/unc, 2. )
          chisqFWHM[j] = chisqFWHM[j]/len(c[0]) 
          chisqBoth[j] = chisqFlux[j] + chisqFWHM[j]
          print "%5.2f  %4.1f  %.2e  %10.2f  %10.2f" % (alpha, r0AU, n0, chisqFlux[j], chisqFWHM[j])
          j = j+1

  # when finished, print out 10 best fits
    ndump = 10
    if ndump > n :
      ndump = n

    nbest = chisqBoth.argsort()[:ndump]
    print "\nbest overall fits:"
    for j in nbest :
      print "%5.2f  %4.1f  %.2e  %10.2f  %10.2f  %10.2f" % (params[j][0], params[j][1], params[j][2], \
        chisqFlux[j], chisqFWHM[j], chisqBoth[j])

    nbest = chisqFlux.argsort()[:ndump]
    print "\nbest Flux fits:"
    for j in nbest :
      print "%5.2f  %4.1f  %.2e  %10.2f  %10.2f  %10.2f" % (params[j][0], params[j][1], params[j][2], \
        chisqFlux[j], chisqFWHM[j], chisqBoth[j])

    nbest = chisqFWHM.argsort()[:ndump]
    print "\nbest FWHM fits:"
    for j in nbest :
      print "%5.2f  %4.1f  %.2e  %10.2f  %10.2f  %10.2f" % (params[j][0], params[j][1], params[j][2], \
        chisqFlux[j], chisqFWHM[j], chisqBoth[j])

  # save results to pickle file for further analysis, if desired
    fout = open("BNmodel.pkl", "wb")
    pickle.dump( [params, chisqFlux, chisqFWHM, chisqBoth], fout )
    fout.close()
    


#  radialModel("BNmodel2.radialModel", 30., -3.5, 7.4, 5.e7, 50, 0.1 )
   

alphaList = numpy.arange(-1.,-4.1,-.1)
r0List = numpy.arange(1.,10.,1.)
n0List = numpy.arange(1.e7,9.e7,1.e7)

def doit () :
  findBestModel( alphaList, r0List, n0List ) 
