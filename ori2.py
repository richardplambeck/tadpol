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
      # allow 1 sec to finish writing tmpfile
    fin = open( tmpfile, "r" )
    for line in fin :
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
    return numpy.split( spectrum, 2, axis=1 )
      # this returns separate freq and flux arrays

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
    fin = open( linelistFile, "r" )
    for line in fin :
      if (len(line) > 1) and (not line.startswith("#")) :
        a = line.split(":")
        print a[2], a[0]
        linefreq = float(a[2])
        species = a[0]
        if (linefreq >= fmin) and (linefreq <= fmax) :
          yval = yvalue( freq, flux, linefreq )
          p.annotate( species, xy=(linefreq,yval), xytext=(linefreq,ymax), horizontalalignment="center", rotation="vertical", size=8, \
             arrowprops=dict( arrowstyle="->", linewidth=0.2 ) )


def hist( freq, flux, plotParams, linelistFile=None ) :
    #pyplot.ioff()
    #pp = PdfPages("spectrum.pdf")

    ymin = plotParams["ymin"]
    ymax = plotParams["ymax"]
    if ymin == ymax :
      ymin,ymax = yminmax( flux )
    print "using ymin,ymax = ",ymin,ymax
   
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
      p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
      fmin = freq[nch1]
      fmax = freq[nch2]
      delta = fmax - fmin
      fmin = fmin - .04*delta
      fmax = fmax + .04*delta
      x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
      p.xaxis.set_major_formatter( x_formatter )
      p.tick_params( axis='both', which='major', labelsize=8 )
      p.axis( [fmin, fmax, ymin, ymax] )
      p.plot( freq[nch1:nch2], flux[nch1:nch2], color="r", linewidth=0.5 )
      ann( p, plotParams["linelistFile"], freq, flux, fmin, fmax, ymax )
      if (np == maxPanelsPerPage) :
        #pyplot.savefig( pp, format='pdf', bbox_inches='tight' )
        pyplot.show()
      nch1 = nch1 + nincr
      nch2 = nch2 + nincr
    if (np != maxPanelsPerPage ) :
      #pyplot.savefig( pp, format='pdf', bbox_inches='tight' )
      pyplot.show()
    #pp.close()

plotParams = { "npanels" : 1,
               "maxPanelsPerPage" : 1,
               "nlap" : 100,
               "ymin" : -.2,
               "ymax" : 2.,
               "linelistFile" : "LineList.csv" }

def doit( infile, region='relpix,box(-2,-2,0,0)', plotParams=plotParams ) :
   [freq, flux] = getspec( infile, region )
   hist( freq, flux, plotParams )
