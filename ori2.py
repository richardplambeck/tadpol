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

clight = 2.99792e5    # speed of light in km/sec

def getspec( infile, region='relpix,box(-2,-2,0,0)', vsource=5., hann=5, tmpfile="junk" ):
    '''dump out spectrum of selected region with imspec, return [chan, freqLSR, flux] arrays'''

  # step 1: use imlist to retrieve velocity and freq information from the header
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

  # step 2: use imspec to dump out the spectrum for the selected region to tmpfile
    chan = []
    freq = []
    flux = []
    p= subprocess.Popen( ( shlex.split("imspec in=%s region=%s options=list,eformat,noheader,hanning,%d log=%s" % \
      (infile,region,hann,tmpfile) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    time.sleep(1)
    result = p.communicate()[0]
    print result
    if "Fatal Error" in result :
      print " --- fatal --- "
      return

  # step 3: read velocities and flux densities from tmpfile, create arrays
    fin = open( tmpfile, "r" )
    for line in fin :
      a = line.split()
      if len(a) > 2 :
        chan.append( int(a[0]) )
        nchan = int( a[0] )
        vlsr = float( a[1] )
        flux.append( float( a[2] ) )
        vlsrcalc = v1 + (nchan - 1) * dv
        if abs(vlsrcalc-vlsr) > 0.05 :
          print "WARNING: channel %d listed vlsr = %.2f, calculated = %.2f" % (nchan,vlsr,vlsrcalc)
        fqLSR = restfreq * (1. - vlsrcalc/clight) 
        freq.append( fqLSR/(1.-vsource/clight) )
        #print nchan, vlsrcalc, fqLSR, freq[-1]
          # freq in rest frame of source
    fin.close()
    print "read in %d lines" % len(freq)

  # step 4: sort in frequency order, return arrays
    spectrum = numpy.array(sorted(zip(freq,chan,flux)))
      # this sorts the chan,freq,flux triplets in freq order
    a,b,c = numpy.split( spectrum, 3, axis=1 )
      # this returns separate freq and flux arrays
    return numpy.reshape(b,len(a)), numpy.reshape(a,len(b)), numpy.reshape(c,len(c))
      # unless I do this the arrays are officially 2D arrays?!

def yminmax( yarray ) :
    '''find ymin and ymax for a plot'''
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

def ann( p, linelistFile, freq, flux, fmin, fmax, ymax ) :
    '''annotates plot with spectral line labels'''
    name,linefreq,TL,intensity,vlsr = processLineList( linelistFile )
    for n in range(0,len(name)) :
        if (linefreq[n] >= fmin) and (linefreq[n] <= fmax) :
          yval = yvalue( freq, flux, linefreq[n] )
          #p.annotate( "%s %.0fK %.3f" % (name[n],TL[n],intensity[n]), xy=(linefreq[n],yval), \
          p.annotate( "%s  %.0fK  %.0fkm/sec" % (name[n],TL[n],vlsr[n]), xy=(linefreq[n],yval), \
             xytext=(linefreq[n],0.95*ymax), horizontalalignment="center", rotation="vertical", size=8, \
             arrowprops=dict( arrowstyle="->", linewidth=0.2 ) )

def findLines( fluxArray, continuum ) :
    '''return convolution of spectrum with lineshape, to locate centers of spectral lines'''
    shape = numpy.zeros(160)
    for n in range( 40,120 ) : 
      shape[n] = 1.
    # shape = signal.gaussian( 100, 10.)
    # return signal.deconvolve( fluxArray-continuum, shape)
    return numpy.convolve( fluxArray-continuum, shape, mode='same')

# plotParams dictionary controls plotting details
plotParams = { "npanels" : 2,
               "maxPanelsPerPage" : 1,
               "nlap" : 100,
               "ymin" : -.5,
               "ymax" : 5.,
               "linelistFile" : "splatalogue.csv", 
               "title" : "Title",      # placeholder for title, to be filled in later
               "contLevel" : 1.5,
               "shadeBad" : False,     # if true, shading shows bad channels, if False shows good channels
               "shadeHgt" : 0.01,       # if shading good channels, shade box ranges from contLevel-shadeHgt to contLevel+shadeHgt
               "convolve" : True,
               "flagtable" : "flags_a0",
               "pdf" : True }

# fluxList is list of 1D amplitude arrays
# each of these is supposed to match the single chan and freq arrays
# the goal is to overlay uvspec or imspec spectra obtained with different resolutions

def hist( chan, freq, fluxList, plotParams, linelistFile=None ) :
    if plotParams["pdf"] :
      pyplot.ioff()
      pp = PdfPages("spectrum.pdf")
    ymin = plotParams["ymin"]
    ymax = plotParams["ymax"]
    if ymin == ymax :
      ymin,ymax = yminmax( fluxList[0] )
      ymax = 1.2*ymax
    print "using ymin,ymax = ",ymin,ymax

  # locate spectral lines, find plot limits
    if plotParams["convolve"] :
      cnv = findLines( fluxList[0], 0.5 )
      cmin,cmax = yminmax( cnv )
      #print "actual cmin,cmax: ",cmin,cmax
      cmin = 1.04 * cmin
      cmax = 1.04 * cmax
      #print "rescaled cmin,cmax: ",cmin,cmax
 
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
      p3 = p.twiny()   
        # p2 shares the same x-axis, p3 the same y-axis
      p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray

      fmin = freq[nch1]
      fmax = freq[nch2]
      delta = fmax - fmin
      fmin = fmin - .04*delta
      fmax = fmax + .04*delta

      chmin = chan[nch1]
      chmax = chan[nch2]
      delta = chmax - chmin
      chmin = chmin - .04*delta
      chmax = chmax + .04*delta

      p.axis( [fmin, fmax, ymin, ymax] )
      if plotParams["convolve"] :
        p2.axis( [fmin, fmax, cmin, cmax] )
        p2.tick_params( axis='y', which='major', labelsize=8, colors='blue' )
        p2.plot( freq[nch1:nch2], cnv[nch1:nch2], color="b", linewidth=0.2, alpha=0.5 )
      p3.axis( [chmin, chmax, ymin, ymax] )

      x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
      p.xaxis.set_major_formatter( x_formatter )
      x_locator = matplotlib.ticker.AutoMinorLocator()

    # this is weird, but the minor axes will appear only on the 2nd of these (whichever their order)
      p3.xaxis.set_minor_locator( x_locator )
      p.xaxis.set_minor_locator( x_locator )

      p.tick_params( axis='both', which='major', labelsize=8 )
      p3.tick_params( axis='x', which='major', labelsize=8 )

      p.plot( freq[nch1:nch2], fluxList[0][nch1:nch2], color="r", linewidth=0.5 )
      if plotParams["contLevel"] :
        p.axhline( y=plotParams["contLevel"], linestyle='--', color='blue')

    # indicate channels that were flagged bad for continuum map
      if plotParams["flagtable"] :
        y0 = numpy.zeros( len(chan) )
        flaggedBad = numpy.zeros( len(chan), dtype=int )
        fin = open( plotParams["flagtable"], "r" )
        for line in fin :
          if not line.startswith("#") :
            a=line.split(",")
            nbstart = int(a[0])
            nbstop = int(a[1])
            for n in range(0, len(chan)) :
              if (chan[n] >= nbstart) and (chan[n] <= nbstop) :
                flaggedBad[n] = 1
        fin.close()
        if plotParams["shadeBad"] :
          p.fill_between( freq, y0, fluxList[0], where=flaggedBad==1, color='red', alpha=0.2 )
        else :
          y1 = plotParams["contLevel"] - plotParams["shadeHgt"]
          y2 = plotParams["contLevel"] + plotParams["shadeHgt"]
          p.fill_between( freq, y1, y2, where=flaggedBad==0, color='blue', alpha=0.1 )

      if plotParams["linelistFile"] :
        ann( p, plotParams["linelistFile"], freq, fluxList[0], fmin, fmax, ymax )
      pyplot.suptitle( plotParams["title"], fontsize=8 )

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

def doit( infile, region='relpix,box(0,-4,2,-2)', vsource=5.0, contLevel=1.7, flagtable="flags_a0", ymin=-.5, ymax=5.,  plotParams=plotParams ) :
# ori2.doit("sp1ch300.cm", region='relpix,box(-172,-5,-170,-3)', flagtable="flags_b1", contLevel=0.11, ymax=0.3, ymin=-.02)  E Src
# ori2.doit("sp2ch300.cm", region='relpix,box(51,-117,53,-115)', flagtable="flags_b2", contLevel=0.065, ymax=0.2, ymin=-.02)  SW Src
# ori2.doit("sp0ch300.cm", flagtable='flags_a0', contLevel=.35, vsource=5.0, region='arcsec,box(4.3,-.2,4.2,-.1)' )
# doit( infile, region='relpix,box(-2,-2,0,0)', vsource=5.0, contLevel=.35, flagtable=None, plotParams=plotParams ) :
# ori2.doit("sp0ch0.cm", flagtable='flags_a0', contLevel=1.55, vsource=5.0, region='relpix,box(-175,-11,-165,-1)')   # E point source
# ori2.doit( "sp0ch250.cm", vsource=7., contLevel=0.05, region="arcsec,box(-1.1,-3.1,-1.4,-2.8)", flagtable="flags_b0")
# def doit( infile, region='relpix,box(151,197,153,199)', plotParams=plotParams ) :
# def doit( infile, region='relpix,box(-30,0,-10,20)', vsource=6., plotParams=plotParams ) :

# def doit( infile, plotParams=plotParams ) :
   [chan, freq, flux ] = getspec( infile, region=region, vsource=vsource )
   plotParams["title"] = infile + " region=" + region + " vsource=%.1d km/sec" % vsource
   plotParams["contLevel"] = contLevel
   plotParams["flagtable"] = flagtable
   plotParams["ymin"] = ymin
   plotParams["ymax"] = ymax
   fluxList = [flux]
   hist( chan, freq, fluxList, plotParams )

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

# read in the lineListFile, in format produced by splatalogue
def processLineList( lineListFile, vsource=5. ) :
    print "reading in spectral line list from %s" % lineListFile
    print "default LSR velocity of spectral lines is %.1f" % vsource
    name = []
    freq = []
    TL = []
    intensity = []
    vdop = []
    dopfac = 1.
    vlsr = vsource
    fin = open( lineListFile, "r" )
    for line in fin :
      if len(line) > 0 :
        a = line.split(":")

      # process header line; figure out what data are in what column
        if (a[0] == "#Species") or (a[0] == "Species") :
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

      # vlsr line gives vlsr of all lines that follow, until next vlsr line
        elif (a[0] == "#VLSR") or (a[0] == "VLSR") : 
          vlsr = float(a[1])
          print "VLSR = %.1f km/sec for following data" % vlsr
          dopfac = 1. - (vlsr - vsource)/clight

      # process data lines
        elif not line.startswith("#") :
          name.append( a[nameCol] )

        # sometimes NRAO-recommended freq corresponds to calculated freq, but other times to measured freq
          if len(a[freqCol]) > 0 :
            freq.append( dopfac * float(a[freqCol]) )
          else :
            freq.append( dopfac * float(a[altFreqCol]) )

          TL.append( float(a[TLcol]) )
          intensity.append( float(a[intensityCol]) )
          vdop.append( vlsr )
          #print "%15s %12.5f %5.0f %6.4f" % ( name[-1], freq[-1], TL[-1], intensity[-1] )  
    fin.close()
    return name,freq,TL,intensity,vdop

def savethese() :
  doit( "sp3ch250.cm", vsource=7., contLevel=0.12, region='arcsec,box(4.46,-.20,4.2,0.04)', flagtable="flags_b3")

