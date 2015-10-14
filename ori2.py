# ori2.py
# non-OO version of ori.py

import math
import time
import cmath
import numpy
import pylab
import sys
import pickle
import subprocess
import shlex
import string
import random
import matplotlib
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from scipy import signal

# dump out spectrum of selected region with imspec, return [freqLSR, flux]
def getspec( infile, region='relpix,box(-2,-2,0,0)', vsource=5., hann=3, tmpfile="junk" ):
    clight = 2.99792e5    # speed of light in km/sec

  # retrieve velocity and freq information from the header
    p= subprocess.Popen( ( shlex.split('imlist in=%s' % infile) ), \
        stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    lines = result.split("\n")
    for line in lines :
      if len(line) > 1 :
        a = line.split()
        n = string.find( line, "restfreq:" )
        if n >= 0 :
          restfreq = float( line[n+9:].split()[0] )
        n = string.find( line, "crval3  :" )
        if n >= 0 :
          v1 = float( line[n+9:].split()[0] )
        n = string.find( line, "cdelt3  :" )
        if n >= 0 :
          dv = float( line[n+9:].split()[0] )
    print "restfreq = %.5f GHz; v1 = %.3f km/sec; dv = %.3f km/sec" % (restfreq,v1,dv)        

  # now dump out the spectrum for the selected region
    freq = []
    flux = []
    p= subprocess.Popen( ( shlex.split("imspec in=%s region=%s options=list,eformat,noheader,hanning,%d log=%s" % \
      (infile,region,hann,tmpfile) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    time.sleep(1)
    result = p.communicate()[0]
      # allow 5 sec to finish writing tmpfile
    fin = open( tmpfile, "r" )
    for line in fin :
      print line
      a = line.split()
      if len(a) > 2 :
        nchan = int( a[0] )
        vlsr = float( a[1] )
        flux.append( float( a[2] ) )
        vlsrcalc = v1 + (nchan-1) * dv
        if abs(vlsrcalc-vlsr) > 0.05 :
          print "WARNING: channel %d listed vlsr = %.2f, calculated = %.2f" % (nchan,vlsr,vlsrcalc)
        fqLSR = restfreq * (1. - vlsrcalc/clight) 
        freq.append( (1.+vsource/clight)*fqLSR )
          # freq in rest frame of source
    fin.close()
    print "read in %d lines" % len(freq)
    spectrum = numpy.array(sorted(zip(freq,flux)))
      # this sorts the freq,flux pairs in freq order
    a,b = numpy.split( spectrum, 2, axis=1 )
      # this returns separate freq and flux arrays
    return numpy.reshape(a,len(a)), numpy.reshape(b,len(b))
      # unless I do this the array is officially a 2D array?!

def yminmax( yarray ) :
    ymin = yarray.min()
    ymax = yarray.max()
    diff = ymax - ymin
    return [ ymin - .05*diff, ymax + .05*diff ]

def yvalue( freq, flux, fline ) :
    '''find amplitude of spectrum at freq fline - use for annotating lines'''
    for n in range(0,len(freq)-1) :
      if (fline >= freq[n]) and (fline < freq[n+1]) : 
        return flux[n]
    return 1000.

# annotates using linelistFile
# linelistFile contains freq, species
def ann( p, linelistFile, freq, flux, fmin, fmax, ymax ) :
    name,linefreq,TL,intensity = processLineList( linelistFile )
    for n in range(0,len(name)) :
        if (linefreq[n] >= fmin) and (linefreq[n] <= fmax) :
          yval = yvalue( freq, flux, linefreq[n] )
          p.annotate( "%s %.0fK %.3f" % (name[n],TL[n],intensity[n]), xy=(linefreq[n],yval), \
             xytext=(linefreq[n],0.95*ymax), horizontalalignment="center", rotation="vertical", size=8, \
             arrowprops=dict( arrowstyle="->", linewidth=0.2 ) )


# slide a template line profile along the spectrum to locate possible spectral lines
# return convolution of freq with lineshape
def findLines( fluxArray, continuum ) :
    shape = numpy.zeros(160)
    for n in range( 40,120 ) : 
      shape[n] = 1.
    # shape = signal.gaussian( 100, 10.)
    # return signal.deconvolve( fluxArray-continuum, shape)
    return numpy.convolve( fluxArray-continuum, shape, mode='same')

def hist( freq, flux, plotParams, linelistFile=None ) :
    if plotParams["pdf"] :
      pyplot.ioff()
      pp = PdfPages("spectrum.pdf")
    ymin = plotParams["ymin"]
    ymax = plotParams["ymax"]
    if ymin == ymax :
      ymin,ymax = yminmax( flux )
      ymax = 1.2*ymax
    print "using ymin,ymax = ",ymin,ymax

  # locate spectral lines, find plot limits
    cnv = findLines( flux, 0.5 )
    cmin,cmax = yminmax( cnv )
    print cmin,cmax
    cmin = -cmax
    cmax = 1.2 * cmax
 
    npanels = plotParams["npanels"]
    maxPanelsPerPage = plotParams["maxPanelsPerPage"]
    nlap = plotParams["nlap"]
    nextra = (npanels+1)*nlap
      # number of frequency channels that are overlaps, plus start and end segments 
    nincr = (len(freq)-nextra)/npanels + nlap
      # nch1 and nch2 increment by this number each time
    nch1 = 0
    nch2 = nincr + nlap - 1
    for npanel in range(0,npanels) :
      np = npanel % maxPanelsPerPage + 1
      p = pyplot.subplot(maxPanelsPerPage,1,np)
      p2 = p.twinx()
        # a twin of Axes sharing the same x-axis
      p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
      fmin = freq[nch1]
      fmax = freq[nch2]
      delta = fmax - fmin
      fmin = fmin - .04*delta
      fmax = fmax + .04*delta
      p.axis( [fmin, fmax, ymin, ymax] )
      p2.axis( [fmin, fmax, cmin, cmax] )
      p2.plot( freq[nch1:nch2], cnv[nch1:nch2], color="b", linewidth=0.2 )
      x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
      p.xaxis.set_major_formatter( x_formatter )
      p.tick_params( axis='both', which='major', labelsize=8 )
      p.plot( freq[nch1:nch2], flux[nch1:nch2], color="r", linewidth=0.5 )
      p.axhline( y=0.52, linestyle='--', color='blue')
      ann( p, plotParams["linelistFile"], freq, flux, fmin, fmax, ymax )
      if (np == maxPanelsPerPage) :
        if plotParams["pdf"] :
          #pyplot.savefig( pp, format='pdf', bbox_inches='tight' )
          pyplot.savefig( pp, format='pdf' )
        pyplot.show()
      nch1 = nch1 + nincr
      nch2 = nch2 + nincr
    if (np != maxPanelsPerPage ) :
      if plotParams["pdf"] :
        # pyplot.savefig( pp, format='pdf', bbox_inches='tight' )
        pyplot.savefig( pp, format='pdf' )
      pyplot.show()
    if plotParams["pdf"] :
      pp.close()

plotParams = { "npanels" : 3,
               "maxPanelsPerPage" : 1,
               "nlap" : 100,
               "ymin" : -.2,
               "ymax" : 2.,
               "linelistFile" : "splatalogue.csv", 
               "pdf" : True }

def doit( infile, region='relpix,box(-2,-2,0,0)', vsource=5.0, plotParams=plotParams ) :
# def doit( infile, region='relpix,box(151,197,153,199)', plotParams=plotParams ) :
# def doit( infile, region='relpix,box(-30,0,-10,20)', vsource=6., plotParams=plotParams ) :
# def doit( infile, plotParams=plotParams ) :
   [freq, flux] = getspec( infile, region=region, vsource=vsource )
   hist( freq, flux, plotParams )

# check definitions of line intensities
# Sij is line strength in Debyes (esu-cm)
# EupperK is energy level of upper level (Kelvin)
# ElowerK is energy level of lower level (Kelvin)
# checked against data for SO+ line at 348115 MHz:
# ... ori2.JPLintensity( 348115., 7.3489, 70.21, 53.51, 624., 300.) gives 0.00773027020636 -2.11180532535
# ...   which seems to match JPL intensity in splatalogue
def JPLintensity( freqMHz, Sba, EupperK, ElowerK, Qrs, T ) :
    # h = 6.626e-27   # cm^2 g/sec
    # c = 2.9979e10   # cm/sec
    # Iba = 8.* pow(math.pi, 3.)/(3.*h*c) * freqMHz * Sba * (math.exp(-EupperK/T) - math.exp(-ElowerK/T))/ Qrs
    Iba = 4.16231e-5 * freqMHz * Sba * (math.exp(-ElowerK/T) - math.exp(-EupperK/T))/ Qrs
    print Iba, math.log10(Iba)

def processLineList( lineListFile ) :
    name = []
    freq = []
    TL = []
    intensity = []
    fin = open( lineListFile, "r" )
    for line in fin :
      if len(line) > 0 :
        a = line.split(":")
        if a[0] == "#Species" :
          nameCol = 0
          for n in range(1,len(a)) :
            if a[n] == "Freq-GHz" :
              freqCol = n
            if a[n] == "E_L (K)" :
              TLcol = n
            if a[n] == "CDMS/JPL Intensity" :
			  intensityCol = n
            if a[n] == "Meas Freq-GHz" :
			  altFreqCol = n
          print "nameCol = %d, freqCol = %d, TLcol = %d, intensityCol = %d" % (nameCol,freqCol, TLcol,intensityCol)
        elif not line.startswith("#") :
          name.append( a[nameCol] )
          if len(a[freqCol]) > 0 :
            freq.append( float(a[freqCol]) )
          else :
            freq.append( float(a[altFreqCol]) )
          TL.append( float(a[TLcol]) )
          intensity.append( float(a[intensityCol]) )
          print "%15s %12.5f %5.0f %6.4f" % ( name[-1], freq[-1], TL[-1], intensity[-1] )  
    fin.close()
    return name,freq,TL,intensity
