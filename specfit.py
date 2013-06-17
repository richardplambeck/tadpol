# specfit.py

import numpy
import re
import scipy.optimize
import mpfit
import mpfitexpr

def gaussian( v, a0, a1, v0, fwhm, slope ) :
  return a0 + a1 * numpy.exp(-2.772 * pow( (v-v0)/fwhm, 2.)) + slope*v 
 
def getdata( infile, xcol, ycol ) :
  x = []
  y = []
  fin = open( infile, "r" ) 
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      x.append( a[xcol-1] )
      y.append( a[ycol-1] )
  fin.close()
  return [numpy.array(x,dtype=float), numpy.array(y,dtype=float)]

def fit( x, y, guess, outfile ) :
  #x,y = getdata( "/fringe2/plambeck/pacs/A86orimsr/H4142narrow.spec", 2, 3 )
  #guess = [ 0.09, 0.047, 23.8, 22.6 ]   
    # continuum level = 0.09
    # line amp = 0.047
    # v0 = 23.8 km/sc
    # FWHM = 22.6 km/sec
  popt,pcov = scipy.optimize.curve_fit( gaussian, x, y, guess )
  print popt
  print pcov
  fout = open( outfile, "w" )
  fout.write("# a0 = %0.4f (%0.4f)\n" % (popt[0],numpy.sqrt(pcov[0][0]) ) )
  fout.write("# a1 = %0.4f (%0.4f)\n" % (popt[1],numpy.sqrt(pcov[1][1]) ) )
  fout.write("# v0 = %0.2f (%0.2f)\n" % (popt[2],numpy.sqrt(pcov[2][2]) ) )
  fout.write("# FWHM = %0.2f (%0.2f)\n" % (popt[3],numpy.sqrt(pcov[3][3]) ) )
  for i in range(0,len(x)) :
    fout.write(" %0.2f  %0.4f  %0.4f\n" % \
	  (x[i],y[i],gaussian(x[i],popt[0],popt[1],popt[2],popt[3]) ) )
  fout.close()
  
def older() :

  x,y = getdata( "/fringe2/plambeck/pacs/A86orimsr/H4142narrow.spec", 2, 3 )
  guess = [ 0.09, 0.047, 23.8, 22.6 ]   
  fit( x, y, guess, "BN_H4142.fit" )

  x,y = getdata( "/fringe2/plambeck/c0676/data/BN.nchUSB.cm.spec", 2, 3 )
  guess = [ 0.24, 0.3, 23.8, 22.6 ]   
  fit( x, y, guess, "BN_H30.fit" )
  
  x,y = getdata( "/fringe2/plambeck/pacs/oripaper/H53alpha.spec", 1, 2 )
  guess = [ 0.24, 0.3, 23.8, 22.6 ]   
  fit( x, y, guess, "BN_H53.fit" )
  
# compute pressure-broadened linewidth in km/sec for H_n_alpha recomb line 
#   if the electron density is ne
# see Brocklehurst and Seaton 1972, eqn 4.6

def pb( n, ne, Te=8000 ) :
  rydberg = 109737.31  # cm-1
  c = 2.998e10  # cm s-1
  delta = 4.7 * pow(n/100., 4.4) * pow(10000./Te, 0.1) * ne  # units: Hz
  freq = rydberg * c * (1./(n*n) - 1./((n+1.)*(n+1)) )
  linewidth = 2.* delta/freq * c/1.e5  # units: km/sec 
  lw2 = c/1.e5 * 1.43e-5 * pow(n/100., 7.4) * pow(10000./Te, 0.1) * (ne/10000.)
  print "freq = %.3f GHz, linewidth = %.1f km/sec" % (freq/1.e9, linewidth)
  print lw2

# ===== end Orion stuff; begin cx355 stuff ===== #

def save() : 
  x,y = getdata( "/fringe2/plambeck/cx355/3c207/src16.dat", 1, 2 )
  guess = [ 0.078, -0.1, 5., 1. ]   
  fit( x, y, guess, "3c207.fit" )
    
# --- this version uses mpfit in order to constrain v0, FWHM --- #

def myfunc( p, fjac=None, x=None, y=None, err=None ) :
  xavg = 0.5 * (x[0] + x[len(x)-1])
  return [0, eval('(y-( p[0] + p[1]*(x-xavg) + p[2]*numpy.exp(-2.772 * pow( (x-p[3])/p[4], 2.)) ) )')]

def newfit() :
  x,y = getdata( "/fringe2/plambeck/cx355/3c207/src16.dat", 1, 2 )
  outfile = '3c207.mfit'
  err = numpy.ones( len(x), dtype=float )
  p0 = [ 0.8,-0.2, 4.5, 1.0, 0. ]   
  parinfo = [ {'fixed':0}, {'fixed':0}, {'fixed':1}, {'fixed':1}, {'fixed':0} ]
    # continuum level = 0.08
    # line amp = -0.02
    # v0 = 5 km/sc
    # FWHM = 1 km/sec
#  [params,yfit] = mpfitexpr( "p[0] + p[1]*numpy.exp(-2.772 * pow( (x-p[2])/p[3], 2.)", x, y, err, p0 ) 
  fa = {'x':x, 'y':y, 'err':err}
  m = mpfit.mpfit(myfunc, p0, functkw=fa, parinfo=parinfo)
  fout = open( outfile, "w" )
  fout.write("# a0 = %0.4f\n" %  m.params[0] )
  fout.write("# a1 = %0.4f\n" %  m.params[1] )
  fout.write("# v0 = %0.2f\n" %  m.params[2] )
  fout.write("# FWHM = %0.2f\n" % m.params[3] )
  for i in range(0,len(x)) :
    fout.write(" %0.2f  %0.4f  %0.4f  %.4f  %.4f\n" % \
	  (x[i],y[i],gaussian(x[i],m.params[0],m.params[1],m.params[2],m.params[3],m.params[4]), \
       y[i] - m.params[4]*(x[i] - m.params[2]) , \
	   gaussian(x[i],m.params[0],m.params[1],m.params[2],m.params[3],0.)+m.params[4]*m.params[2]) ) 
  fout.close()
 
# fit linear slope + gaussian to spectrum
# spectrum is in columns xcol and ycol of infile
# x-units can be chan or vlsr
# y-units can be Jy, K, or absorption depth
# parameter list: 
#   yavg = avg continuum level at center of spectrum (in y-units)
#   yslope = slope (in y-units/x-units)
#   amp = amplitude of gaussian (in y-units)
#   v0 = center of gaussian (in x-units)
#   fwhm = fwhm of gaussian (in x-units)
# inputs:
#   guess = list of starting values [avg, slope, amp, v0, fwhm]
#   constrain = list of dictionary values specifying whether values are fixed or not
#      (default constraint is to let all values vary)
# output file:
#   col 1 = x
#   col 2 = y
#   col 3 = yfit 
#   col 4 = y with linear slope removed
#   col 5 = yfit with linear slope removed (i.e., Gaussian only)

def newfit2( infile, outfile, xcol, ycol, p0, \
    constrain=[ {'fixed':0}, {'fixed':0}, {'fixed':0}, {'fixed':0}, {'fixed':0} ] ) :
  [x,y] = getdata( infile, xcol, ycol )
  err = numpy.ones( len(x), dtype=float )
  parinfo = [ {'fixed':0}, {'fixed':0}, {'fixed':0}, {'fixed':0}, {'fixed':0} ]
  fa = {'x':x, 'y':y, 'err':err}
  m = mpfit.mpfit(myfunc, p0, functkw=fa, parinfo=constrain)

  fout = open( outfile, "w" )
  fout.write("# %s\n" %  outfile )
  fout.write("# fit to spectrum %s, xcol=%d, ycol=%d\n" % ( infile, xcol, ycol ) )
  fout.write("# final fit:\n")
  fout.write("#   yavg= %0.5f   %d \n" %  (m.params[0], constrain[0]['fixed'])  )
  fout.write("#   slope = %0.5f   %d\n" % (m.params[1], constrain[1]['fixed'])  )
  fout.write("#   amp = %0.5f   %d\n" %  (m.params[2], constrain[2]['fixed'])  )
  fout.write("#   v0 = %0.3f   %d\n" %  (m.params[3], constrain[3]['fixed'])  )
  fout.write("#   FWHM = %0.3f   %d\n" % (m.params[4], constrain[4]['fixed'])  )
  fout.write("#\n")

  xavg = 0.5 * (x[0] + x[len(x)-1])
  fit = eval('m.params[0] + m.params[1]*(x-xavg) + \
     m.params[2]*numpy.exp(-2.772 * pow( (x-m.params[3])/m.params[4], 2.))')
     # fit to raw data, including slope
  yminusSlope = eval('y - m.params[1]*(x-xavg)')
     # raw data with linear slope removed
  gfit = eval('m.params[0] + m.params[2]*numpy.exp(-2.772 * pow( (x-m.params[3])/m.params[4], 2.))')
     # gaussian fit to data with slope removed
  yminusSlopeNorm = eval('yminusSlope/m.params[0]')
     # normalized absorption 
  gfitNorm = eval('gfit/m.params[0]')
  for i in range(0,len(x)) :
    fout.write(" %10.3f  %9.4f  %9.4f  %9.4f  %9.4f %9.5f %9.5f\n" % \
	  (x[i],y[i],fit[i],yminusSlope[i],gfit[i],yminusSlopeNorm[i],gfitNorm[i] ) )
  fout.close()
 
# boxcar average Bon-Chul's spectra, and convert to Jy
def binit( infile, outfile, nbin, Jy_source ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  xavg = yavg = 0.
  navg = 0
  for line in fin :
    a = line.split()
    xavg = xavg + float(a[0])
    yavg = yavg + float(a[1])
    navg = navg + 1
    if navg >= nbin :
      xavg = xavg/navg
      yavg = yavg/navg
      fout.write("%11.4f %11.4f %11.5f\n" % ( xavg, yavg, Jy_source + 15.*yavg ) )
      xavg = 0.
      yavg = 0.
      navg = 0
  fin.close()
  fout.close()

def doit() :

  p0 = [ 5.5, 0., -.2, -9.7, 1.0 ]   
  newfit2( "3C454.3_dickplot.txt", "3c454_OH.spec", 2, 3, p0)
 
  v0 = -9.54      # based on fit to OH
  fwhm = 1.26     # based on fit to OH
  p0 = [ 5.5, 0., -.2, v0, fwhm ]   
  constrain = [ {'fixed':0}, {'fixed':0}, {'fixed':0}, {'fixed':1}, {'fixed':1} ]

  newfit2( "/fringe2/plambeck/cx355/3c454/spectra/c16.dat", "3c454_CO.spec", 1, 4, p0, constrain=constrain)

  newfit2( "/fringe2/plambeck/cx355/3c454/spectra/c15.dat", "3c454_CN.spec", 1, 4, p0, constrain=constrain)

  newfit2( "/fringe2/plambeck/cx355/3c454/spectra/c8.dat", "3c454_CS.spec", 1, 4, p0, constrain=constrain)

  binit( "HCO+Abs_3C454.3_YS.dat", "tmp", 3, 5.5 )
  newfit2( "tmp", "3c454_HCO+.spec", 1, 3, p0, constrain=constrain)

  binit( "HCNAbs_3C454.3_YS.dat", "tmp", 8, 5.5 )
  newfit2( "tmp", "3c454_HCN.spec", 1, 3, p0, constrain=constrain)

