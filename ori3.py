# ori2.py
# non-OO version of ori.py

import math
import time
import cmath
import numpy
import sys
import pickle
import subprocess
import shlex
import string
import random
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import signal

ckms = 2.99792e5    # speed of light in km/sec


# ----------------------------------------------------------------------------------------------
# Note: - axis3 of spectral line data cube is LSR velocity relative to restfreq (at LSR=0)
#       - from this, one can compute the restfreq of that channel in the LSR=0 frame (the lab frame) 
#       - getspec divides by doppler factor to compute restfreq in the source frame (lsr=vsource) 
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

      # just to be safe, compare vlsr from file with vlsr computed from crval3
        vlsrcalc = v1 + (nchan - 1) * dv
        if abs(vlsrcalc-vlsr) > 0.05 :
          print "WARNING: channel %d listed vlsr = %.2f, calculated = %.2f" % (nchan,vlsr,vlsrcalc)

        fqLSR = restfreq * (1. - vlsrcalc/ckms) 
          # ... a signal at restfreq emitted by gas at vlsrcalc would be observed at fqLSR;
          #     in other words, this is the lab freq in the LSR=0 frame
        freq.append( fqLSR/(1.-vsource/ckms) )
          # ... line frequency in the rest frame of the source!

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


# getuvspec is retrieves visibility spectra using uvspec

def getuvspec( infile, vsource=5., hann=3, tmpfile="junk", select='uvrange(0,200)' ):
    '''dump out spectrum of selected region with uvspec, return [chan, freqLSR, flux] arrays'''

  # step 1: use uvlist to retrieve veldop from the header
    veldop = 0.
    p= subprocess.Popen( ( shlex.split('uvlist vis=%s options=var,full' % infile) ), \
        stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    lines = result.split("\n")
    for line in lines :
      if line.startswith("veldop") :
        a = line.split(":")
        veldop = float( a[1] )
    print "veldop = %.3f km/sec" % veldop       

  # step 2: use uvspec to dump out the spectrum for the selected region to tmpfile
    chan = []
    freq = []
    flux = []
    p= subprocess.Popen( ( shlex.split("uvspec vis=%s select=%s axis=freq,amp stokes=i \
      options=ampscalar,avall,nobase interval=100000 hann=%d log=%s" % \
      (infile,select,hann,tmpfile) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    time.sleep(1)
    result = p.communicate()[0]
    print result
    if "Fatal Error" in result :
      print " --- fatal --- "
      return

  # step 3: read velocities and flux densities from tmpfile, create arrays
    fin = open( tmpfile, "r" )
    nline = 0
    for line in fin :
      nline = nline + 1
      a = line.split()
      chan.append( nline )
      freq.append( float( a[0] )/(1.-(vsource-veldop)/ckms) )
        # correct to doppler frame of source
      flux.append( float( a[1] ) )
    fin.close()
    print "read in %d lines" % len(freq)

  # step 4: sort in frequency order, return arrays
    spectrum = numpy.array(sorted(zip(freq,chan,flux)))
      # this sorts the chan,freq,flux triplets in freq order
    a,b,c = numpy.split( spectrum, 3, axis=1 )
      # this returns separate freq and flux arrays
    return numpy.reshape(b,len(a)), numpy.reshape(a,len(b)), numpy.reshape(c,len(c))
      # unless I do this the arrays are officially 2D arrays?!

# ----------------------------------------------------------------------------------------------
# processLineList reads in file produced by splatalogue
# the goal is to compare lab frequencies from splatalogue with observed frequencies that
#   are reported in the vlsr=vsource frame
# so if line is emitted by gas at LSR velocity vsource, return lab freq directly
# if line is emitted at a different LSR velocity, doppler correct the frequency
#
# used by two routines:
#   - ann uses it to annotate the plot (all outputs needed)
#   - pruneLineList uses it to read in the raw file (variables veldop and shadevel not used)

def processLineList( lineListFile, vsource=5. ) :
    '''read splatalogue.csv file, return [ name, freq, TL, intensity, veldop, shadevel ] lists'''

    print "reading in spectral line list from %s" % lineListFile
    name = []
    freq = []
    QN = []		# note: this routine does not return QN currently
    TL = []
    intensity = []
    vdop = []
    shadeVel = []

    vlsr = vsource   
    print "default LSR velocity of spectral lines is %.1f" % vsource
      # ... the default vlsr is vsource, in which case no doppler correction will be applied

    fin = open( lineListFile, "r" )
    for line in fin :
      if len(line) > 0 :
        a = line.split(":")

      # process header line; figure out what data are in what column
        if "Species" in a[0] :
          nameCol = 0
          for n in range(1,len(a)) :
            if (a[n] == "Freq-GHz") or ("Ordered" in a[n]) :
              freqCol = n
            if "Meas Freq-GHz" in a[n] :
			  altFreqCol = n
            if "E_L (K)" in a[n] :
              TLcol = n
            if "CDMS/JPL Intensity" in a[n] :
			  intensityCol = n
            if "Resolved QNs" in a[n] :
			  QNcol = n
          print "nameCol = %d, freqCol = %d, TLcol = %d, QNcol = %d, intensityCol = %d" % \
             (nameCol,freqCol,QNcol,TLcol,intensityCol)

      # vlsr line gives vlsr of all lines that follow, until next vlsr line
        elif "VLSR" in line : 
          vlsr = float(a[1])
          print "VLSR = %.1f km/sec for data that follow" % vlsr

      # process data lines
        elif not line.startswith("#") :
          name.append( a[nameCol] )
          vdop.append( vlsr )
          dopfac = 1. - (vlsr - vsource)/ckms
            # ... doppler correct lab freq if line is emitted at velocity different than vsource

        # sometimes freq shows up as calculated freq, other times as measured freq
          if len(a[freqCol]) > 0 :
            freq.append( dopfac * float(a[freqCol]) )
          else :
            freq.append( dopfac * float(a[altFreqCol]) )

          QN.append( a[QNcol] )
          TL.append( float(a[TLcol]) )
          intensity.append( float(a[intensityCol]) )
          #print "%15s %12.5f %5.0f %6.4f" % ( name[-1], freq[-1], TL[-1], intensity[-1] )  
 
        # lines to be shaded are indicated by "SHADE", then vlsr to be shaded in next column
          shade = 0.
            # ... default is not to shade the line
          for i in range( 0, len(a)-1 ) :
            if "SHADE" in a[i] :
              shade = float( a[i+1] )
          shadeVel.append( shade )

    fin.close()
    return name,freq,TL,intensity,vdop,shadeVel

# ----------------------------------------------------------------------------------------------
# parse splatalogue output, remove unobserved frequencies, duplicates

def pruneLineList( infile, outfile ) :
    freqRanges = [ [340.545, 344.310], [352.665, 356.430], [649.282, 651.177], [661.303, 665.067], [665.922, 667.817] ]
    copy = []
    name = []
    QN = []
    freq = []
    cat = []

  # pass 1: read file, create lists of species, freq, QN, linelist for each input line 
    fin = open( infile, "r" )
    for line in fin :
      copy.append( 1 )     # default is to copy this line
      name.append( "" )      
      freq.append( 0. )
      QN.append( "" )
      cat.append( "" )
      if len(line) > 0 :
        a = line.split(":")
        if "Species" in a[0] :
          nameCol = 0
          for n in range(1,len(a)) :
            if "Ordered" in a[n] :
              freqCol = n
            if "Resolved" in a[n] :
              QNcol = n
            if "Linelist" in a[n] :
              catCol = n
        else : 
          name[-1] = a[nameCol]
          freq[-1] = float(a[freqCol])
          QN[-1] = a[QNcol]
          cat[-1] = a[catCol]

        # for catalog lines, copy only those in observed freq ranges
          copy[-1] = 0
          for fr in freqRanges :
            if (freq[-1] >= fr[0]) and (freq[-1] <= fr[1]) :
              copy[-1] = 1           
    fin.close()

  # now weed out duplicate lines; favor CDMS > JPL > SLAIM
    for n in range(0,len(copy)) :
      if copy[n] and name[n] :         # don't fiddle with header, or blank, or comment lines
        nsrchmax = n+50				# no point in searching more than 50 lines ahead
        if nsrchmax > len(copy) :
          nsrchmax = len(copy)
        for m in range(n+1,nsrchmax) :
          if (name[m] == name[n] ) and (QN[m].replace(" ","") == QN[n].replace(" ","")) :
            print freq[m],cat[m]
            if "CDMS" in cat[n] :
              copy[m] = 0
            elif "CDMS" in cat[m] :
              copy[n] = 0 
            elif "SLAIM" in cat[n] :
              copy[n] = 0 
            elif "SLAIM" in cat[m] :
              copy[m] = 0
            
  # pass 2: open file again, write out valid lines to output file
    fin = open( infile, "r" )
    fout = open( outfile, "a" )
    n = 0
    for line in fin :
      if copy[n] :
        fout.write("%s:SHADE:40.\n" % line.rstrip('\n'))
      n = n+1
    fin.close()
    fout.close()
      
# ----------------------------------------------------------------------------------------------
def yminmax( yarray ) :
    '''find ymin and ymax for a plot'''
    ymin = yarray.min()
    ymax = yarray.max()
    diff = ymax - ymin
    return [ ymin - .05*diff, ymax + .05*diff ]

# ----------------------------------------------------------------------------------------------
def yvalue( freq, flux, fline ) :
    '''find amplitude of spectrum at freq fline - use for annotating lines'''
    for n in range(0,len(freq)-1) :
      if (fline >= freq[n]) and (fline < freq[n+1]) : 
        return flux[n]
    return 1000.

# ----------------------------------------------------------------------------------------------
def ann( p, freq, flux, fmin, fmax, plotParams ) :
    '''annotates plot with spectral line labels'''
    if plotParams["linelistFile"] :
      name,linefreq,TL,intensity,vlsr,shadeVel = processLineList( plotParams["linelistFile"] )
      for n in range(0,len(name)) :
        if (linefreq[n] >= fmin) and (linefreq[n] <= fmax) :
          yval = yvalue( freq, flux, linefreq[n] )
          p.annotate( "%s  %.0fK  %.0fkm/sec" % (name[n],TL[n],vlsr[n]), xy=(linefreq[n],yval), \
             xytext=(linefreq[n],0.95*plotParams["ymax"]), horizontalalignment="center", rotation="vertical", size=8, \
             arrowprops=dict( arrowstyle="->", linewidth=0.2 ) )
          if shadeVel[n] != 0. :
            shadeColor = "yellow"
            if shadeVel[n] < 0. :
              shadeColor = "blue"
            deltaGHz = abs(shadeVel[n])/(2.*2.998e5) * linefreq[n]
            p.fill_between( freq, plotParams["contLevel"], flux, color=shadeColor, alpha=1.0, \
              where=numpy.abs(freq-linefreq[n]) < deltaGHz )
            
             
# ----------------------------------------------------------------------------------------------
#    contLevel,rms = findRms( chan, flux, flagtable )
def findRMS( chan, flux, flagtable ) :
    flaggedBad = numpy.zeros( len(chan), dtype=int )
    if flagtable :
      try :
        fin = open( flagtable, "r" )
      except :
        print "couldn't find flagtable %s - skipping rms and contLevel calculation" % flagtable
        return [0.,0.,flaggedBad] 

      for line in fin :
        if not line.startswith("#") :
          a=line.split(",")
          nbstart = int(a[0])
          nbstop = int(a[1])
          for n in range(0, len(chan)) :
            if (chan[n] >= nbstart) and (chan[n] <= nbstop) :
              flaggedBad[n] = 1
      fin.close()
      maskedFlux = numpy.ma.array( flux, mask=flaggedBad )
      contLevel = numpy.ma.mean( maskedFlux )
      rms = numpy.ma.std( maskedFlux )
      print "calculated contLevel = %.2f" % contLevel
      print "calculated rms = %.3f" % rms
      return [contLevel,rms,flaggedBad] 
    

# ----------------------------------------------------------------------------------------------
def findLines( fluxArray, contLevel, rms, chansToConvolve ) :
    '''return convolution of spectrum with lineshape, to locate centers of spectral lines'''

  # normFlux is excess over continuum expressed in terms of rms
    normFlux = (fluxArray-contLevel)/(math.sqrt(chansToConvolve)*rms)

    shape = numpy.ones(chansToConvolve)
    # shape = signal.gaussian( 100, 10.)
    # return signal.deconvolve( fluxArray-continuum, shape)
    return numpy.convolve( normFlux, shape, mode='same')

# plotParams dictionary controls plotting details
plotParams = { "npanels" : 10,
               "maxPanelsPerPage" : 1,
               "nlap" : 100,
               "ymin" : -.5,
               "ymax" : 5.,
               #"linelistFile" : "/o/plambeck/OriALMA/Spectra/splatalogue.csv", 
               "linelistFile" : "/o/plambeck/Downloads/splatalogue.csv", 
               "title" : "Title",      # placeholder for title, to be filled in later
               "contLevel" : 0.,
               "rms" : 0.0,       # if shading good channels, shade box ranges from contLevel-shadeHgt to contLevel+shadeHgt
               "convolve" : True,
               "flagtable" : None,
               "pdf" : True }

# ----------------------------------------------------------------------------------------------
# fluxList is list of 1D amplitude arrays
# each of these is supposed to match the single chan and freq arrays
# the goal is to overlay uvspec or imspec spectra obtained with different resolutions

def hist( chan, freq, fluxList, flaggedBad, plotParams ) :

    if plotParams["pdf"] :
      pyplot.ioff()
      pp = PdfPages("spectrum.pdf")
    ymin = plotParams["ymin"]
    ymax = plotParams["ymax"]
    if ymin == ymax :
      ymin,ymax = yminmax( fluxList[0] )
      ymax = 1.2*ymax
    plotParams["ymin"] = ymin
    plotParams["ymax"] = ymax

  # locate spectral lines, find plot limits
    if plotParams["convolve"] :
      chansToConvolve = round( (36./3.e5)*freq[ len(freq)/2 ]/(freq[1]-freq[0]) )
         # expect emission lines from SrcI to be 36 km/sec wide
      print "convolve %d chans" % chansToConvolve
      cnv = findLines( fluxList[0], plotParams["contLevel"], plotParams["rms"], chansToConvolve )
      cmin,cmax = yminmax( cnv )
      cmin = 1.04 * cmin
      cmax = 1.04 * cmax
 
  # default is to plot 350 GHz spectra in 2 panels because of large number of spectral lines
    npanels = plotParams["npanels"]
    if npanels == 0 :
      if freq[0] < 400. : npanels = 2
      else : npanels = 1

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

    # indicate channels that were NOT flagged bad for continuum map
      if plotParams["flagtable"] :
        y1 = plotParams["contLevel"] - plotParams["rms"]
        y2 = plotParams["contLevel"] + plotParams["rms"]
        p.fill_between( freq, y1, y2, where=flaggedBad==0, color='0.9' )

    # annotate lines in linelistFile; also, if desired, shade emission or absorption
      if plotParams["linelistFile"] :
        ann( p, freq, fluxList[0], fmin, fmax, plotParams )

      pyplot.suptitle( plotParams["title"], fontsize=10 )
      if (np == maxPanelsPerPage) :
        if plotParams["pdf"] :
          pyplot.savefig( pp, format='pdf' )
        pyplot.show()
      nch1 = nch1 + nincr
      nch2 = nch2 + nincr
    if (np != maxPanelsPerPage ) :
      if plotParams["pdf"] :
        pyplot.savefig( pp, format='pdf' )
      pyplot.show() 
    if plotParams["pdf"] :
      pp.close()

# ori2.doit("sp1ch300.cm", region='relpix,box(-172,-5,-170,-3)', flagtable="flags_b1", contLevel=0.11, ymax=0.3, ymin=-.02)  E Src
# ori2.doit("sp2ch300.cm", region='relpix,box(51,-117,53,-115)', flagtable="flags_b2", contLevel=0.065, ymax=0.2, ymin=-.02)  SW Src
# ori2.doit("sp0ch300.cm", flagtable='flags_a0', contLevel=.35, vsource=5.0, region='arcsec,box(4.3,-.2,4.2,-.1)' )
# doit( infile, region='relpix,box(-2,-2,0,0)', vsource=5.0, contLevel=.35, flagtable=None, plotParams=plotParams ) :
# ori2.doit("sp0ch0.cm", flagtable='flags_a0', contLevel=1.55, vsource=5.0, region='relpix,box(-175,-11,-165,-1)')   # E point source
# ori2.doit( "sp0ch250.cm", vsource=7., contLevel=0.05, region="arcsec,box(-1.1,-3.1,-1.4,-2.8)", flagtable="flags_b0")
# def doit( infile, region='relpix,box(151,197,153,199)', plotParams=plotParams ) :
# def doit( infile, region='relpix,box(-30,0,-10,20)', vsource=6., plotParams=plotParams ) :

      
# ----------------------------------------------------------------------------------------------
def doit( infile, region='relpix,box(0,-4,2,-2)', vsource=5.0, contLevel=1.7, flagtable="flags_a0", ymin=-.5, ymax=5.,  plotParams=plotParams ) :

  # fill in plotParams from input keywords 
    plotParams["title"] = infile + " region=" + region + " vsource=%.1d km/sec" % vsource
    plotParams["contLevel"] = contLevel
    plotParams["flagtable"] = flagtable
    plotParams["ymin"] = ymin
    plotParams["ymax"] = ymax

  # retrieve the spectrum
    if infile.endswith( ".uv" ) :
      [chan, freq, flux ] = getuvspec( infile, vsource=vsource )
    else :
      [chan, freq, flux ] = getspec( infile, region=region, vsource=vsource )

  # compute contLevel and rms from unflagged channels
    if flagtable :
      [contLevel, rms, flaggedBad] = findRMS( chan, flux, flagtable )
      plotParams["contLevel"] = contLevel
      plotParams["rms"] = rms
    else :
      plotParams["contLevel"] = 0.
      plotParams["rms"] = 0. 
      flaggedBad = []

    fluxList = [flux]
    hist( chan, freq, fluxList, flaggedBad, plotParams )


# ----------------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------------
def savethese() :
  doit( "sp3ch250.cm", vsource=7., contLevel=0.12, region='arcsec,box(4.46,-.20,4.2,0.04)', flagtable="flags_b3")

# plot visibilities (raw, model, or residuals) in 3D
# infile must be logfile created by uvlist
def visplot( infile ) :
  uin = []
  vin = []
  ampin = []
  header = True
  fin = open(infile, "r")
  for line in fin:
    if header :
      if "Vis" in line : 
        header = False
    else :
      a = line.split()                   
      if float(a[7]) > 0.0 :
        uin.append( float(a[5]) )
        vin.append( float(a[6]) )
        ampin.append( float(a[7]) )
  fin.close()
  u = numpy.array(uin)
  v = numpy.array(vin)
  amp = numpy.array(ampin)
  fig = pyplot.figure()
  ax = Axes3D(fig)
  ax.scatter(u,v,amp,c=amp,cmap=cm.rainbow)
  #fig.colorbar(q)
  pyplot.show()

# this is a helper file for checkresid, below
# it reads log created by uvlist, returns real and imaginary parts of each visibility
def readvis( infile ) :
  ampin = []
  phsin = []
  header = True
  fin = open(infile, "r")
  for line in fin:
    if header :
      if "Vis" in line : 
        header = False
    else :
      a = line.split()                      # split line into string tokens
      ampin.append( float(a[7]) )
      phsin.append( math.pi/180.*float(a[8]) )
  fin.close()
  amp = numpy.array(ampin)
  phs = numpy.array(phsin)
  return [amp*numpy.cos(phsin), amp*numpy.sin(phsin)]

# I had doubts about whether uvfit returned correct residuals, so wrote this to check
# this is a very fragile routine, depends on the 3 input files containing exactly the
#    same number of lines, with listings in exactly the same order... but it was
#    adequate for a quick check
def checkresid( vis, model, resid, outfile ) :
  vr,vi = readvis( vis )
  mr,mi = readvis( model )
  rr,ri = readvis( resid )
  fout = open( outfile, "w")
  for n in range(0,len(vr)) :
    fout.write("%8.4f %8.4f  %8.4f %8.4f\n" % (vr[n]-mr[n], vi[n]-mi[n], rr[n], ri[n]))
  fout.close()

# returns line intensities in LINEAR units for a particular temperature T
def calcIntensity( infile, T ) :
  T0 = 300.
  hPlanck = 6.62607e-27 
  kBoltzmann = 1.38065e-16
  name,freq,TL,intensity,vdop,shadeVel = processLineList( infile )
  Iba = numpy.zeros( len(name) )
  for n in range(0, len(name) ) :
    TU = TL[n] + (hPlanck * freq[n]*1.e9 / kBoltzmann)   # upper state energy in K
    Iba[n] = pow(10.,intensity[n]) * \
      (math.exp(-TU/T) - math.exp(-TL[n]/T)) / (math.exp(-TU/T0) - math.exp(-TL[n]/T0))
  return Iba

def plotIntensity( infile ) :
  name,freq,TL,intensity,vdop,shadeVel = processLineList( infile )
  T = 100
  Iba = calcIntensity( infile, T) * 37424./pow(T,1.5)
  fig,ax = pyplot.subplots()
  pyplot.bar( freq, Iba, .01, color="red", edgecolor="red" )
  pyplot.show()

# snippet extracts sections of spectra centered on lines in splatalogue list
# create a 'Line' dictionary for each one
# pickle the 'Line' dictionaries one by one, append to outfile

def snippet( infile, outfile, region='relpix,box(-1,-1,1,1)', vsource=5.0, flagtable="flags_a0", \
     linelistFile="/o/plambeck/OriALMA/Spectra/splatalogue.csv", vrange=80. ) :

  # retrieve the spectrum
    if infile.endswith( ".uv" ) :
      [chan, freq, flux ] = getuvspec( infile, vsource=vsource )
    else :
      [chan, freq, flux ] = getspec( infile, region=region, vsource=vsource )

  # compute contLevel and rms from unflagged channels
    if flagtable :
      [contLevel, rms, flaggedBad] = findRMS( chan, flux, flagtable )
    else :
      contLevel = 0.

    name,linefreq,TL,intensity,vlsr,shadeVel = processLineList( linelistFile )
    for n in range(0,len(name)) :
      FillingIn = False
      Line = {}
      for i in range(0,len(freq)-1) :
        dv = (freq[i]-linefreq[n])/linefreq[n] * ckms   # velocity relative to line center (if line is at nominal VLSR)
        if abs(dv) < vrange/2. :
          if not FillingIn :
            FillingIn = True
            Line["name"] = name[n]
            Line["TL"] = TL[n]
            Line["linefreq"] = linefreq[n]
            Line["intensity"] = intensity[n]
            Line["contLevel"] = contLevel
            vlsr0 = vlsr[n]
            Line["vel"] = []
            Line["flux"] = []
          Line["vel"].append( dv + vlsr0 )
          Line["flux"].append( flux[i] )
        else :     # past end of snippet
          if FillingIn :
            print "dumping to pickle", Line
            fout = open( outfile, "ab" )     # binary because I'll be pickling to it; append because I'll append
            pickle.dump( Line, fout )
            fout.close() 
            FillingIn = False
      if FillingIn :     # past end of spectrum
        print "dumping to pickle", Line
        fout = open( outfile, "ab" )     # binary because I'll be pickling to it; append because I'll append
        pickle.dump( Line, fout )
        FillingIn = False
        fout.close() 

# plot the snippets 

def plotSnippets( infile ) :
  npanels = 1
  maxPanelsPerPage = 1
  Lines = []   # list of Line dictionaries
  fin = open( infile, "rb" )
  while (True) :
    try :
      Line = pickle.load( fin )
      Lines.append( Line )
    except :
      break 
  print "read in %d Line objects" % ( len(Lines) )

  for npanel in range(0,npanels) :
    np = npanel % maxPanelsPerPage + 1
    p = pyplot.subplot(maxPanelsPerPage,1,np)
    p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
    vmin = -60.
    vmax = 70.
    p.axis( [-60., 70., -0.2, 8.] )
    off = 0.
    for Line in Lines :
      offset = off * numpy.ones(len(Line["flux"]))
      p.plot( Line["vel"], Line["flux"]+offset, color="r", linewidth=1. )
      p.axhline( y=Line["contLevel"]+off, linestyle='--', color='blue')
      off = off + .5
    pyplot.show()

# special-purpose routine to read fluxes from xxGHz.csv files and plot
csvFileList = [ '90GHz.csv', '229GHz.csv', '350GHz.csv', '660GHz.csv' ]
srcNameList = [ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'BN' ]
def plotFluxes( ) :
  nplot=0
  for srcName in srcNameList :
    nplot = nplot+1
    p = pyplot.subplot(3,5,nplot)
    p.axis( [70., 1000., 2., 10000.] )
    p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
    frqGHz = []
    SmJy = []
    uncSmJy = []
    print " "
    for csvFile in csvFileList :
      fin = open( csvFile, 'r' )
      for line in fin :
        if len(line) > 0 :
          a = line.split(',')
          if (len(a[0]) > 0 ) and (a[0] == srcName) :
            try :
              mJy = 1000. * float(a[10]) * float(a[16]) 
              if mJy > 0. :
                n = csvFile.find("GHz") 
                frqGHz.append( float( csvFile[0:n] ) )
                SmJy.append( mJy )
                unc = 1000. * float(a[11]) * float(a[16]) 
                uncSmJy.append( math.sqrt( pow(unc,2.) + pow(0.1*SmJy[-1],2.) ) )
              else :
                pass     # don't include if flux is zero
            except :             # handles case where flux = '-'
              pass
    if len(frqGHz) == len(SmJy) :
      Sf2 = []
      for frq,S in zip(frqGHz,SmJy) : 
        if frq > 200. and frq < 240. :
          fref = frq
          Sref = S
      for frq,S in zip(frqGHz,SmJy) : 
          Sf2.append( pow(frq/fref, 2.) * Sref )
      p.loglog( frqGHz, SmJy )
      p.loglog( frqGHz, Sf2, linestyle='dashed')
      p.errorbar( frqGHz, SmJy, yerr=uncSmJy)
      p.text( .1, .8, "%s" % srcName, transform=p.transAxes )
      p.tick_params( axis='both', which='major', labelsize=8 )
    else :
      print srcName, frqGHz, SmJy
  pyplot.show()


