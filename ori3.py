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
import os
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import signal

ckms = 2.99792e5    # speed of light in km/sec

# plotParams dictionary controls plotting details
plotParams = { "npanels" : 2,
               "maxPanelsPerPage" : 2,
               "nlap" : 20,
               "ymin" : -.2,       # -.1 for band7, -.2 for band9
               "ymax" : 5.8,      # 1.95 for band7, 5.5 for band9
               "linelistFile" : "/o/plambeck/OriALMA/Spectra/splat_ann.csv", 
               "title" : "Title",      # placeholder for title, to be filled in later
               "contLevel" : 0.,
               "rms" : 0.0,       # if shading good channels, shade box ranges from contLevel-shadeHgt to contLevel+shadeHgt
               "convolve" : False,
               "flagtable" : None,
               "shadeUnflagged" : False,	  # show which channels were unflagged?
               "labelChanAxis" : False,        # show channel axis along top of plots?
               "pdf" : True }

# create a "specList" dictionary for use by snippet2 and pubFig
#   'file' is text file containing the spectrum integrated over some region
#   'flag' is flags file indicating which spectra were flagged for continuum maps
B7spw0 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw0.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b0"}
B7spw1 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw1.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b1"}
B7spw2 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw2.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b2"}
B7spw3 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw3.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b3"}
B9spw0 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw0.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a0"}
B9spw1 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw1.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a1"}
B9spw2 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw2.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a2"}
B9spw3 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw3.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a3"}
specList = [B7spw0,B7spw1,B7spw2,B7spw3,B9spw0,B9spw1,B9spw2,B9spw3]


# run getspec, dump spectrum to a text file; adjust freq scale to be correct for source at LSR velocity vsource
def dumpspec( infile, outfile, region='arcsec,box(.1)', vsource=0., hann=3 ):
    [chan, freq, flux ] = getspec( infile, region=region, vsource=vsource )
    fout = open( outfile, "w" )
    fout.write("# infile: %s\n" % infile)
    fout.write("# region: %s\n" % region)
    fout.write("# vsource: %.2f\n" % vsource)
    fout.write("# hann: %d\n" % hann)
    fout.write("# chan, freq, flux\n")
    for [chn,frq,flx] in zip( chan, freq, flux ) :
      fout.write("%5d  %12.6f  %9.7f\n" % (chn,frq,flx))
    fout.close() 

# retrieve spectrum from a text file; adjust freq scale to be correct for source at LSR velocity vsource
def readspec( inspec, vsource=5. ) :
    chan = []
    freq = []
    flux = []
    fin = open( inspec, "r")
    for line in fin :
      a = line.split()
      if "vsource:" in line :
        vs0 = float(a[2])
      if not line.startswith("#") :
        chan.append( int(a[0]) )
        freq.append( float(a[1])/(1.-(vsource-vs0)/ckms) )
        flux.append( float(a[2]) ) 
    fin.close()
    return numpy.array(chan), numpy.array(freq), numpy.array(flux)

# ----------------------------------------------------------------------------------------------
# getspec returns spectrum of selected region from spectral line cube; note that region
#  command could also select channel range, if desired
# Note: - axis3 of spectral line data cube is LSR velocity relative to restfreq (at LSR=0)
#       - from this, one can compute the restfreq of that channel in the LSR=0 frame (the lab frame) 
#       - getspec divides by doppler factor to compute restfreq in the source frame (lsr=vsource) 

def getspec( infile, region='relpix,box(-2,-2,0,0)', vsource=0., hann=3, tmpfile="junk" ):
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
      # this returns separate freq and flux array
    return numpy.reshape(b,len(a)), numpy.reshape(a,len(b)), numpy.reshape(c,len(c))
      # unless I do this the arrays are officially 2D arrays?!


# ----------------------------------------------------------------------------------------------
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
# processLineList reads in file in splatalogue.csv format
# the goal is to compare lab frequencies from splatalogue with observed frequencies that
#   are reported in the vlsr=vsource frame
# so if line is emitted by gas at LSR velocity vsource, return lab freq directly; this
#   is the usual situation
# it is possible to add "VLSR:xx:" lines to splatalogue.csv; if such a line is read, all
#    subsequent data (until another VLSR line is read) will be Doppler corrected by
#    1 - (vlsr-vsource)/c frame; IN THIS CASE THE FREQUENCY WILL NOT BE THE LAB FREQUENCY
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
    TU = []     # not returned 
    Smusq = []  # not returned
    intensity = []
    vdop = []
    shadeColor = []     # color to shade line
    shadeWidth = []     # shade this width in km/sec
    colorCol = widthCol = None

    vlsr = vsource   
    print "default LSR velocity of spectral lines is %.1f" % vsource
      # ... the default vlsr is vsource, in which case no doppler correction will be applied

    fin = open( lineListFile, "r" )
    #freqCol = altFreqCol = TLcol = TUcol = intensityCol = SmusqCol = QNcol = None
    for line in fin :
      if len(line) > 1 :
        a = line.split(":")

      # process header line; figure out what data are in what column
        if "Species" in a[0] :
          SmusqCol = TLcol = TUcol = intensityCol = None
          nameCol = 0
          for n in range(1,len(a)) :
            if (a[n] == "Freq-GHz") or ("Ordered" in a[n]) :
              freqCol = n
            if "Meas Freq-GHz" in a[n] :
			  altFreqCol = n
            if "E_L (K)" in a[n] :
              TLcol = n
            if "E_U (K)" in a[n] :
              TUcol = n
            if "CDMS/JPL Intensity" in a[n] :
			  intensityCol = n
            if "S<sub>ij</sub>&#956;" in a[n] :
              SmusqCol = n
            if "Resolved QNs" in a[n] :
			  QNcol = n
            if "Width" in a[n] :
              widthCol = n
            if "Color" in a[n] :
              colorCol = n
          #print "nameCol = %d, freqCol = %d, TLcol = %d, QNcol = %d, intensityCol = %d, colorCol = %d, widthCol = %d" % \
          #   (nameCol,freqCol,QNcol,TLcol,intensityCol,colorCol,widthCol)

      # vlsr line gives vlsr of all lines that follow, until next vlsr line
        elif "VLSR" in line : 
          vlsr = float(a[1])
          print "VLSR = %.1f km/sec for data that follow" % vlsr

      # process data lines; this can crash if one of the required columns is not present
        elif not line.startswith("#") :
          print line
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

          if SmusqCol :
            Smusq.append( float(a[SmusqCol]) )
          else :
            Smusq.append( 0. )

          if TUcol :
            TU.append( float(a[TUcol]) )
          else :
            TU.append( -1. )

          if TLcol :
            TL.append( float(a[TLcol]) )
          else :
            TL.append( -1. )

          if intensityCol :
            intensity.append( float(a[intensityCol]) )
          else :
            intensity.append( 0. )
          

        # colorCol and widthCol will be None if we are processing raw splatalogue file
          if colorCol :
            shadeColor.append( a[colorCol] )
          else :
            shadeColor.append( None )
          if widthCol :
            shadeWidth.append( float(a[widthCol]) )
          else :
            shadeWidth.append( None )

    fin.close()
    return name,freq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth

# ----------------------------------------------------------------------------------------------
# parse splatalogue output, remove unobserved frequencies, duplicates 

def pruneLineList( infile="/o/plambeck/Downloads/splatalogue.csv", outfile="/o/plambeck/OriALMA/Spectra/splat_ann.csv" ) :
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
            
  # pass 2: open file again, append valid lines to output file
    fin = open( infile, "r" )
    fout = open( outfile, "a" )
    n = 0
    for line in fin :
      if copy[n] :
        fout.write("%s:None:None:\n" % line.rstrip('\n'))
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
      if ((fline >= freq[n]) and (fline < freq[n+1])) or ((fline <= freq[n]) and (fline > freq[n+1])) : 
        return flux[n]
    return 1000.

# ----------------------------------------------------------------------------------------------
def ann( p, freq, flux, fmin, fmax, plotParams ) :
    '''annotates plot with spectral line labels'''
    print "annotating freq range %.3f to %.3f GHz" % (fmin,fmax)
    if plotParams["linelistFile"] :
      name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( plotParams["linelistFile"] )
      for n in range(0,len(name)) :
        if (linefreq[n] >= fmin) and (linefreq[n] <= fmax) :
          print "fmin = %.4f, fmax = %.4f, linefreq = %.4f" % (fmin,fmax,linefreq[n])
          yval = yvalue( freq, flux, linefreq[n] )
          if yval < plotParams["contLevel"] : yval = plotParams["contLevel"]
          if yval > 0.8*plotParams["ymax"] : yval = plotParams["contLevel"]
          if TU[n] > 0. :
            p.annotate( "%s  %.0fK" % (name[n],TU[n]), xy=(linefreq[n],yval), \
             xytext=(linefreq[n],0.92*plotParams["ymax"]), horizontalalignment="center", rotation="vertical", size=7, \
             arrowprops=dict( arrowstyle="-", linewidth=0.2 ) )    # use arrowstyle="->" to get arrow
          else :
            p.annotate( "%s" % (name[n]), xy=(linefreq[n],yval), \
             xytext=(linefreq[n],0.92*plotParams["ymax"]), horizontalalignment="center", rotation="vertical", size=7, \
             arrowprops=dict( arrowstyle="-", linewidth=0.2 ) )
          if shadeColor[n] :
            deltaGHz = shadeWidth[n]/(2.*2.998e5) * linefreq[n]
            fillmax = linefreq[n]+deltaGHz
            if fillmax > fmax : fillmax = fmax
            fillmin = linefreq[n]-deltaGHz
            if fillmin < fmin : fillmin = fmin
            midpt = (fillmax + fillmin)/2.
            diff = (fillmax - fillmin)/2.
            alpha = 0.7 
            if shadeColor[n] == "green" : alpha=0.5
            if shadeColor[n] == "blue" : alpha=0.2
            p.fill_between( freq, plotParams["contLevel"], flux, color=shadeColor[n], alpha=alpha, lw=0, \
              where=numpy.abs(freq-midpt) < diff )
            
             
# ----------------------------------------------------------------------------------------------
# return continuum level, continuum rms, and flaggedBad array
# DANGER - spectra are expected to cover the entire spectral window; partial spectra will mess up

def findRMS( chan, freq, flux, flagtable ) :

    flaggedBad = numpy.zeros( len(chan), dtype=int )
    if flagtable :
      try :
        fin = open( flagtable, "r" )
      except :
        print "couldn't find flagtable %s - skipping rms and contLevel calculation" % flagtable
        return [0.,0.,flaggedBad] 

    # note that flagtable is based on 3840 chans/spectrum, whereas spectrum we are
    #   evaluating may have been made with line=chan,3840/nscale,1,nscale,nscale
    # figure out nscale from number of channels in spectrum; then, to prevent mistakes,
    #   confirm that frequency increment makes sense
    # note that there is no protection against mixing up flagtables
      
      nscale = 3840/len(chan)
      freqStep = (freq[-1] - freq[0])/len(freq)
      expectedFreqStep = nscale * 1.875/3840.
      if abs(freqStep - expectedFreqStep) > 1.e-5 :
        print "ERROR: incorrect freqStep - skippin grms and contLevel calculation"
        return [0.,0.,flaggedBad]
 
      for line in fin :
        if not line.startswith("#") :
          a=line.split(",")
          nbstart = int(a[0])
          nbstop = int(a[1])
          for n in range(0, len(chan)) :
            if (nscale*chan[n] >= nbstart) and (nscale*chan[n] <= nbstop) :
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

# ---------------------------------------------------------
# draw dotted boxes around good channel ranges
# flaggedBad is array; good channels are 0, bad are 1

def boxBase( ax, y1, y2, freq, flaggedBad ) :
  Good = False    
  goodStart = 0
  goodStop = 3839
    # ... Good = True indicates that we are processing a Good range (part of future rectange)
  for n in range(0,len(flaggedBad)) :
    print "n = %d, flaggedBad = %d, Good = %d, GoodStart = %d, GoodStop = %d" % (n, flaggedBad[n], Good, goodStart, goodStop)
    if Good and flaggedBad[n] :   # new bad chan range is starting; finish Good rectangle, if any
      goodStop = n-1
      goodRect = Rectangle( (freq[goodStart],y1), (freq[goodStop]-freq[goodStart]), (y2-y1), \
        fill=True, edgecolor='red', facecolor='yellow', alpha=1.)   # target range
      ax.add_patch( goodRect )
      Good = False
      print "new good range: %d to %d" % (goodStart,goodStop)
    elif not Good and not flaggedBad[n] :   # new good chan range is starting
      goodStart = n
      Good = True
  if Good :
    goodStop = n
    goodRect = Rectangle( (freq[goodStart],y1), (freq[goodStop]-freq[goodStart]), (y2-y1), \
        fill=True, edgecolor='red', facecolor='yellow', alpha=0.8, linewidth=1 )   # target range
    ax.add_patch( goodRect )
 
# ----------------------------------------------------------------------------------------------
# fluxList is list of 1D amplitude arrays
# each of these is supposed to match the single chan and freq arrays
# the goal is to overlay uvspec or imspec spectra obtained with different resolutions

def hist( chan, freq, fluxList, flaggedBad, plotParams ) :
    colorTable = ["red","blue","orange","black"]
    if plotParams["pdf"] :
      pyplot.ioff()
      pp = PdfPages("spectrum.pdf")
    else :
      pyplot.ion()
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
      # p2 = p.twinx()
        # p2 shares the same x-axis, p3 the same y-axis
      p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
      if npanel == npanels-1 :
         p.set_xlabel("freq (GHz)", fontsize=8)
      p.set_ylabel("flux density (Jy)", fontsize=8)

      fmin = freq[nch1]
      print "nch1 = %d, fmin = %.5f" % (nch1,fmin)
      fmax = freq[nch2]
      delta = fmax - fmin
      fmin = fmin - .04*delta
      fmax = fmax + .04*delta

      chmin = chan[nch1]
      chmax = chan[nch2]
      delta = chmax - chmin    # note that delta can be negative!
      chmin = chmin - .04*delta
      chmax = chmax + .04*delta

      p.axis( [fmin, fmax, ymin, ymax] )
      # if plotParams["convolve"] :
      #   p2.axis( [fmin, fmax, cmin, cmax] )
      #   p2.tick_params( axis='y', which='major', labelsize=8, colors='blue' )
      #   p2.plot( freq[nch1:nch2], cnv[nch1:nch2], color="b", linewidth=0.2, alpha=0.5 )
      if plotParams["labelChanAxis"] :
        p3 = p.twiny()   
        p3.axis( [chmin, chmax, ymin, ymax] )
        p3.xaxis.set_minor_locator( x_locator )
        p3.tick_params( axis='x', which='major', labelsize=8 )

      x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
      p.xaxis.set_major_formatter( x_formatter )
      x_locator = matplotlib.ticker.AutoMinorLocator()

    # this is weird, but the minor axes will appear only on the 2nd of these (whichever their order)
      p.xaxis.set_minor_locator( x_locator )

      p.tick_params( axis='both', which='major', labelsize=8 )

      for flux,color in zip(fluxList,colorTable) :
        p.plot( freq[nch1:nch2], flux[nch1:nch2], color=color, linewidth=0.5 )
      if plotParams["contLevel"] :
        p.axhline( y=plotParams["contLevel"], linestyle='--', color='blue')

    # indicate channels that were NOT flagged bad for continuum map
      if plotParams["flagtable"] and plotParams["shadeUnflagged"] :
        y1 = plotParams["contLevel"] - plotParams["rms"]
        y2 = plotParams["contLevel"] + plotParams["rms"]
        p.fill_between( freq, y1, y2, where=flaggedBad==0, color='0.9' )

      # alternative for Fig 1 of paper: draw dotted rectanges around good ranges
      #  boxBase( p, y1, y2, freq, flaggedBad )      

    # annotate lines in linelistFile; also, if desired, shade emission or absorption
      if plotParams["linelistFile"] :
        ann( p, freq, fluxList[0], freq[nch1], freq[nch2], plotParams )

      pyplot.suptitle( plotParams["title"], fontsize=8 )
      if (np == maxPanelsPerPage) :
        if plotParams["pdf"] :
          pyplot.savefig( pp, format='pdf' )
        pyplot.show()
      nch1 = nch1 + nincr
      nch2 = nch2 + nincr

  # finish up; plot any leftover panels
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
    [chan, freq, flux ] = getspec( infile, region=region, vsource=vsource )

  # compute contLevel and rms from unflagged channels
    [contLevel, rms, flaggedBad] = findRMS( chan, freq, flux, flagtable )
    plotParams["contLevel"] = contLevel
    plotParams["rms"] = rms

    fluxList = [flux]
    hist( chan, freq, fluxList, flaggedBad, plotParams )

# -----------------------------------------------------------------------------------------------
# this is a faster version that retrieves spectrum from text file, rather than running imspec each time
def doit2( inspec, inspec2=None, inspec3=None, vsource=5.0, contLevelOverride=0.,\
   flagtable="flags_b1", ymin=-.1, ymax=3.,  plotParams=plotParams ) :

  # fill in plotParams from input keywords 
    #plotParams["title"] = inspec + "(red), " + inspec2 + "(blu), " + inspec3 + "(org), vsource=%.1d km/sec" % vsource
    plotParams["title"] = inspec + "  vsource=%.1d km/sec" % vsource
    plotParams["flagtable"] = flagtable
    plotParams["ymin"] = ymin
    plotParams["ymax"] = ymax

  # retrieve the spectrum
    [chan, freq, flux ] = readspec( inspec, vsource=vsource )
    if (inspec2) :
      [chan2, freq2, flux2 ] = readspec( inspec2, vsource=vsource )
    if (inspec3) :
      [chan3, freq3, flux3 ] = readspec( inspec3, vsource=vsource )

  # compute contLevel and rms from unflagged channels
    [contLevel, rms, flaggedBad] = findRMS( chan, freq, flux, flagtable )
    plotParams["contLevel"] = contLevel
    plotParams["rms"] = rms
    
    fluxList = [flux]
    if inspec2 : fluxList.append( flux2 )
    if inspec3 : fluxList.append( flux3 )
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

# designed to make figures for publication
# each figure is a full page with 4 panels, each covering a single spectral window

def pubFig( specList=[B9spw0,B9spw1,B9spw2,B9spw3], plotParams=plotParams ) :
#def pubFig( specList=[B7spw0,B7spw1,B7spw2,B7spw3], plotParams=plotParams ) :

    vsource=5.
    if plotParams["pdf"] :
      pyplot.ioff()
      pp = PdfPages("spectrum.pdf")
    ymin = plotParams["ymin"]
    ymax = plotParams["ymax"]
    npanels = 4
    fig = pyplot.figure( figsize=(8,11) )

    npanel = 1
    for spectrum in specList: 
      [chan, freq, flux ] = readspec( spectrum["file"], vsource=vsource )
      p = pyplot.subplot(npanels,1,npanel)
      p.grid( True, linewidth=0.05, color="0.01" )   # color=0.1 is a light gray

      fmin = freq[0]
      fmax = freq[-1]
      delta = fmax - fmin
      fmin = fmin - .04*delta
      fmax = fmax + .04*delta

      p.axis( [fmin, fmax, ymin, ymax] )
      x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
      p.xaxis.set_major_formatter( x_formatter )
      x_locator = matplotlib.ticker.AutoMinorLocator()
      p.xaxis.set_minor_locator( x_locator )
      p.tick_params( axis='both', which='major', labelsize=8 )
      p.plot( freq, flux, color="r", linewidth=0.5 )
      [contLevel, rms, flaggedBad] = findRMS( chan, freq, flux, spectrum["flag"] )
      plotParams["contLevel"] = contLevel
      p.axhline( y=plotParams["contLevel"], linestyle='--', color='blue')

    # label the spw
      ii = string.find( spectrum["file"], "spw" )
      p.text( 1.03, .91, "%s" % spectrum["file"][ii:ii+4], transform=p.transAxes, \
        horizontalalignment='right', fontsize=7, rotation='vertical' )

    # annotate lines in linelistFile; also, if desired, shade emission or absorption
      if plotParams["linelistFile"] :
        ann( p, freq, flux, freq[0], freq[-1], plotParams )
      if npanel == npanels :
         p.set_xlabel("freq (GHz)", fontsize=8)
      p.set_ylabel("flux density (Jy)", fontsize=8)
      npanel = npanel + 1

    if plotParams["pdf"] :
      pyplot.savefig( pp, format='pdf' )
      pp.close()
    pyplot.show()

# returns line intensities in LINEAR units for a particular temperature T
def calcIntensity( infile, T ) :
  T0 = 300.
  hPlanck = 6.62607e-27 
  kBoltzmann = 1.38065e-16

  name,freq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( infile )
  Iba = numpy.zeros( len(name) )
  for n in range(0, len(name) ) :
    TU = TL[n] + (hPlanck * freq[n]*1.e9 / kBoltzmann)   # upper state energy in K
    Iba[n] = pow(10.,intensity[n]) * \
      (math.exp(-TU/T) - math.exp(-TL[n]/T)) / (math.exp(-TU/T0) - math.exp(-TL[n]/T0))
  return Iba

def plotIntensity( infile ) :
  name,freq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( infile )
  T = 100
  Iba = calcIntensity( infile, T) * 37424./pow(T,1.5)
  fig,ax = pyplot.subplots()
  pyplot.bar( freq, Iba, .01, color="red", edgecolor="red" )
  pyplot.show()

# snippet extracts sections of spectra centered on lines listed in "snip.csv" (splatalogue format)
# create a 'Line' dictionary for each one
# pickle the 'Line' dictionaries one by one, append to outfile
# note that desired velocity range could extend off the end of the spectrum, in which case the
#   integrated fluxes will be lower limits; for now, check for this rare circumstance manually
# although code allows one to extract snippet from map or even from raw visibility data, I
#   am not actively using this feature - so be very careful if using these alternative methods
# I am assuming that spectra in the txt files were created using Hanning=3, so that number
#   of independent channels is n/2; thus estimated uncertainty in red and blue flux densities
#   is rms/sqrt(nred/2) or rms/sqrt(nblue/2), where rms is computed from the unflagged channels
#   of the entire spectrum

def snippet2( lineListFile="/o/plambeck/OriALMA/Spectra/snip.csv",
    specList=specList, outfile="snip",
    vrange=65., vsource=5.0 ) :

  # read the lineList file
    name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( lineListFile )

  # extract snippet for each line in the file
    for n in range(0,len(name)) :
 
    # find spectral window that contains the line
      for spectrum in specList :

      # retrieve the spectrum from txt file, from uv data, or from map (txt is preferred)
        if spectrum["file"].endswith( ".txt" ) :
          [chan, freq, flux ] = readspec( spectrum["file"], vsource=vsource )
        elif spectrum["file"].endswith( ".uv") :
          [chan, freq, flux ] = getuvspec( infile, vsource=vsource )
        else :
          [chan, freq, flux ] = getspec( infile, region=region, vsource=vsource )

      # extract the snippet if the line is in this spectral window
        if (linefreq[n] > freq[0]) and (linefreq[n] < freq[-1]) :
          print "found %.3f GHz in file %s" % (linefreq[n],spectrum["file"])
          [contLevel, rms, flaggedBad] = findRMS( chan, freq, flux, spectrum["flag"] )
          FillingIn = False
          Line = {}
          for i in range(0,len(freq)-1) :
            dv = (1. - freq[i]/linefreq[n]) * ckms   # velocity relative to line center (if line is at nominal VLSR)
            if abs(dv) < vrange/2. :
              if not FillingIn :
                FillingIn = True
                Line["name"] = name[n]
                Line["QN"] = QN[n]
                Line["TL"] = TL[n]
                Line["TU"] = TU[n]
                Line["Smusq"] = Smusq[n]
                Line["linefreq"] = linefreq[n]
                Line["intensity"] = intensity[n]
                Line["contLevel"] = contLevel
                Line["shadeColor"] = shadeColor[n]
                Line["shadeWidth"] = shadeWidth[n]
                vlsr0 = vdop[n]
                Line["vel"] = []
                Line["flux"] = []
                Line["freq"] = []
                Line["intfluxRed"] = 0.
                Line["intfluxBlue"] = 0.
                nred = 0
                nblue = 0
              vel = dv + vlsr0
              Line["vel"].append( vel )
              Line["flux"].append( flux[i] )
              Line["freq"].append( freq[i] )
              if (vel > 5.) and (vel <= 23.) :
                Line["intfluxRed"] = Line["intfluxRed"] + (flux[i] - contLevel)
                nred = nred + 1
                # print "red wing :  %.3f  %.3f" % (vel, flux[i]-Line["contLevel"] )
              if (vel > -13.) and (vel <= 5.) :
                Line["intfluxBlue"] = Line["intfluxBlue"] + (flux[i] - contLevel)
                nblue = nblue + 1
                # print "blue wing:  %.3f  %.3f" % (vel, flux[i]-Line["contLevel"] )

            else :     # past end of snippet
              if FillingIn :
                Line["rmsRed"] = rms * math.sqrt(nred/2.)     # assume that spectra were Hanning smoothed! 
                Line["rmsBlue"] = rms * math.sqrt(nblue/2.)   # assume that spectra were Hanning smoothed! 
                print "rmsRed,Blue = ",Line["rmsRed"],Line["rmsBlue"]
                fout = open( outfile, "ab" )     # binary because I'll be pickling to it; append because I'll append
                pickle.dump( Line, fout )
                fout.close() 
                FillingIn = False
          if FillingIn :     # past end of spectrum
            Line["rmsRed"] = rms * math.sqrt(nred/2.)     # assume that spectra were Hanning smoothed! 
            Line["rmsBlue"] = rms * math.sqrt(nblue/2.)   # assume that spectra were Hanning smoothed! 
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
      off = off + 3.
    pyplot.show()

# special-purpose routine to read fluxes from xxGHz.csv files and plot
# csvFileList = [ '90GHz.csv', '229GHz.csv', '350GHz.csv', '660GHz.csv' ]
csvFileList = [ '229GHz_500.csv', '350GHz_500.csv', '660GHz_500.csv' ]
srcNameList = [ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' ]
def plotFluxes( ) :
  nplot=0
  for srcName in srcNameList :
    nplot = nplot+1
    p = pyplot.subplot(2,4,nplot)
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

def makeEllintList( csvFile ) :
  fin = open( csvFile, 'r' )
  n = csvFile.find("GHz") 
  frqGHz = csvFile[0:n] 
  for line in fin :
    if len(line) > 0 :
      a = line.split(',')
      if (len(a) > 16 ) :
        print "# source %s" % a[0]
        print "ellint in=$FILE center=%s,%s radius=.025,1.,.025 scale=%s log=ellint_%s_%s" % \
          (a[14],a[15],a[16],frqGHz,a[0])
  fin.close()

def getuvamps( infile, select, tmpfile="junk" ):
    '''dump out visibility amplitudes'''
    u = []
    v = []
    amp = []
    if select == None :
      select = "uvrange(0,10000)"
    p= subprocess.Popen( ( shlex.split("uvlist vis=%s select=%s recnum=1000000 log=%s" % \
      (infile,select,tmpfile) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    time.sleep(1)
    result = p.communicate()[0]
    print result
    if "Fatal Error" in result :
      print " --- fatal --- "
      return
    fin = open( tmpfile, "r" )
    for line in fin :
      a = line.split()
      if (len(a) > 8) and ("XX" in a[4] or "YY" in a[4] or "I" in a[4] or "LL" in a[4] or "RR" in a[4] ) :
        u.append( float(a[5]) )
        v.append( float(a[6]) )
        amp.append( float(a[7]) )
    fin.close()
    return numpy.array(u), numpy.array(v), numpy.array(amp)

# plot visibilities (raw, model, or residuals) in 3D
# infile is a visibility file with XX, YY (won't work for LL, RR)
def visplot( infile, select=None ) :
  u,v,amp = getuvamps( infile, select=select )
  fig = pyplot.figure()
  ax = Axes3D(fig)
  ax.scatter(u,v,amp,c=amp,cmap=cm.rainbow)
  ax.scatter(-u,-v,amp,c=amp,cmap=cm.rainbow)
  ax.auto_scale_xyz( [-900,900], [-900,900], [0,2] )
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

def visrms( infile, select=None ) :
  u,v,amp = getuvamps( infile, select )
  print numpy.std(amp)

def TbPixel( Jy, sizeArcsec, freqGHz ) :
  clight = 2.998e10
  hPlanck = 6.62e-27
  kB = 1.38e-16
  dOmega = pow( sizeArcsec*math.pi/(180.*60.*60.), 2.)
  B = Jy*1.e-23/dOmega   # brightness in cgs units
  exponent = math.log(1. + 2.*hPlanck*pow(freqGHz*1.e9,3.)/(B * pow(clight,2.)))
  Tb = hPlanck*freqGHz*1.e9/(kB * exponent)
  print "Tb = %.1f" % Tb
  return Tb

def TbDisk( Jy, majArcsec, minArcsec, freqGHz ) :
  area = math.pi * majArcsec * minArcsec / 4.
  sizeArcsec = math.sqrt(area)
  Tb = TbPixel( Jy, sizeArcsec, freqGHz)
  actualDiam = majArcsec * 410. * 1.5e13
     # disk diameter in cm
  pradiated = 2. * math.pi * pow(actualDiam/2.,2.) * 5.67e-5 * pow(Tb,4.) / 2.e33
     # power radiated from top and bottom of disk in solar luminosities
  print "luminosity = %.2e Lo" % pradiated
  
# processes file produced by casaplotms
# averages XX and YY; could hanning smooth or resort, if desired
def casaAmp( infile, outfile, vlsr=5 ) :
  frq = numpy.zeros( 3840, dtype=float)
  amp = numpy.zeros( 3840, dtype=float)
  npts = numpy.zeros( 3840, dtype=int)
  fin = open( infile, "r" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      nchn = int(a[2])
      if (nchn < 0) or (nchn > 3839) :
        print "error: nchn = %d", nchn
      else :
        frq[nchn] = frq[nchn] + float(a[10])
        amp[nchn] = amp[nchn] + float(a[1])
        npts[nchn] = npts[nchn]+1
  fin.close()
  print "npts:", npts
  fout = open( outfile, "w" )
  fout.write("# output of ori3.casaAmp( '%s', '%s', vlsr=%.2f\n" % (infile,outfile,vlsr) )
  fout.write("# vlsr is the assumed source velocity in km/sec\n")
  fout.write("# chan  freq(VLSR=0.00) freq(VLSR=%.2f)    amp\n" % vlsr)
  for n in range(0, 3840) :
    f = frq[n]/npts[n]
    fvlsr = f * (1.+ vlsr/2.998e5)
       # this is the rest freq of the line that would be observed in this 
       # channel if it were emitted by a source moving at VLSR=vlsr km/sec
    fout.write("%4d  %14.9f  %14.9f %10.6f\n" % (n, f, fvlsr, amp[n]/npts[n]) )
  fout.close()

# reads miriad style flag file with channels numbered from 1-3840;
# subtracts 1 to put in casa style with chans numbered 0-3839, then appends into big list
# NOTE: BUG - last semicolon in each list needs to be replaced by comma
def casaFlags( flagFileList ) :
  n = 0
  fout = open("casaFlags","w")
  for flagFile in flagFileList :
    fout.write("%d:" % n)
    n = n+1
    fin = open( flagFile, "r" )
    for line in fin :
      a = line.split(',')
      fout.write("%d~%d;" % (int(a[0])-1,int(a[1])-1) )
    fin.close()
  fout.close()
  
# dump out text files with spectra
# "a" series is an 0.2 arcsec box
# "b" series is an 0.3 arcsec box

regiona = "arcsec,box(.1)"
regionb = "arcsec,box(.15)"
regionc = "arcsec,box(-.3,-.45,-.45,-.3)"
regiond = "arcsec,box(1.3,2,1,2.3)"

def makespec() :
  dumpspec( "spw0.100plus.cm", "spw0.100plus.c.txt", region=regionc, vsource=0., hann=3 )
  dumpspec( "spw1.100plus.cm", "spw1.100plus.c.txt", region=regionc, vsource=0., hann=3 )
  dumpspec( "spw2.100plus.cm", "spw2.100plus.c.txt", region=regionc, vsource=0., hann=3 )
  dumpspec( "spw3.100plus.cm", "spw3.100plus.c.txt", region=regionc, vsource=0., hann=3 )

# helper routine for evalSnippets; returns math.log(value), or -100 if value is negative
def evalLog( value ) :
  if value > 0. :
    return math.log(value)
  else :
    return -100.

# make table of integrated intensities for rotational diagram using data in 'snip' file
def rotdiag( infile="snip", outfile="rotdiag" ) :
    Lines = []
    fin = open( infile, "rb" )
    while (True) :
      try :
        Line = pickle.load( fin )
        Lines.append( Line )
      except :
        break 
    fin.close()
    print "read in %d Line objects" % ( len(Lines) )
   
    fout = open( outfile, "w" ) 
    fout.write("# ori3.rotdiag; using red+blue wings\n")
    for Line in Lines :
      Sdv = Line["intfluxRed"] + Line["intfluxBlue"]
      unc = 1.4*(Line["rmsBlue"])   # very approximate!
      #Sdv = Line["intfluxBlue"]
      #unc = Line["rmsBlue"]
      c = 1./(Line["linefreq"]*Line["Smusq"])
      print " %11.5f   %.4e %.4e   %.4e %.4e %.4e  %8.1f" % \
       ( Line["linefreq"], Sdv, unc, evalLog(c*Sdv), evalLog(c*(Sdv-unc)), evalLog(c*(Sdv+unc)), Line["TU"] )
      fout.write("# %11.5f   %.4e %.4e   %.4e %.4e %.4e  %8.1f\n" % \
       ( Line["linefreq"], Sdv, unc, evalLog(c*Sdv), evalLog(c*(Sdv-unc)), evalLog(c*(Sdv+unc)), Line["TU"] ) )
      fout.write("move %.2f %.4e\n" % (Line["TU"], evalLog(c*Sdv) ) )
      fout.write("dot\n")
      fout.write("move %.2f %.4e\n" % (Line["TU"], evalLog(c*(Sdv-unc)) ) )
      fout.write("draw %.2f %.4e\n" % (Line["TU"], evalLog(c*(Sdv+unc)) ) )
    fout.close()
     
# read rotdiag.csv file, open .mc maps to get mean brightness temp in 0.2" box, write wip commands
def rotdiag2( inlistFile="/o/plambeck/OriALMA/Spectra/rotdiag.csv", outfile="rotdiagGuts.wip", vrange=[5.,20.] ) :

    kB = 1.38e-16    # Boltzman's constant in cgs units

  # read list of lines that must be processed
    name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( inlistFile )

  # open wip file, write out header info
    fout = open( outfile, "w" ) 
    fout.write("# created using ori3.rotdiag2(); using velocity range %.1f to %.1f\n" % (vrange[0],vrange[1]))

  # process lines one at a time; find "mc" map, compute integrated intensity
    for n in range(0,len(linefreq)) :
      mcFile = findpvFile( linefreq[n], searchType=".mc" ) 
      if mcFile :
        Sdv,nch = getW( mcFile, vrange )   # get integrated emission; units must be K-km/sec

      # uncertainty estimates are based on BLANK_354.25.mc (~1 K rms) and BLANK_649.80.mc (~ 3K rms);
      #   assume that uncertainty in integrated Tb is sqrt(nch) * unc/chan
        rms = 1.
        if linefreq[n] > 500. : rms = 3.
        unc = math.sqrt(nch) * rms
        print "W = %.2f unc = %.2f" % (Sdv,unc)
        const = 3.* kB * 1.e5/(8. * pow(math.pi,3) * linefreq[n] * 1.e9 * Smusq[n] * pow(1.e-18, 2) )
          # formula A4 in Cummins et al 1986; 1.e5 is to convert km/sec to cm/sec; 1.e-36 converts debye^2 to (esu-cm)^2
        print " %11.5f   %.4e %.4e   %.4e %.4e %.4e  %8.1f" % \
         ( linefreq[n], Sdv, unc, evalLog(const*Sdv), evalLog(const*(Sdv-unc)), evalLog(const*(Sdv+unc)), TU[n] )
        fout.write("# %11.5f   %.4e %.4e   %.4e %.4e %.4e  %8.1f\n" % \
         ( linefreq[n], Sdv, unc, evalLog(const*Sdv), evalLog(const*(Sdv-unc)), evalLog(const*(Sdv+unc)), TU[n]) )
        fout.write("move %.2f %.4e\n" % (TU[n], evalLog(const*Sdv) ) )
        fout.write("dot\n")
        fout.write("move %.2f %.4e\n" % (TU[n], evalLog(const*(Sdv-unc)) ) )
        fout.write("draw %.2f %.4e\n" % (TU[n], evalLog(const*(Sdv+unc)) ) )
      else :
        print "skipping line at freq %.3f - couldn't find file %s" % (linefreq[n],mcFile)
    fout.close()

# helper routine that opens "mc" map file (line minus continuum), uses imspec to compute mean Tb in region
def getW( mcFile, vrange, region="arcsec,box(0.1)", tmpfile="junk" ) :

    p= subprocess.Popen( ( shlex.split("imspec in=%s region=%s options=list,noheader,tb plot=mean log=%s" % \
      (mcFile,region,tmpfile) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    time.sleep(1)
    result = p.communicate()[0]
    print result   # to spot errors

  # read velocities and temperatures from tmpfile, compute integrted intensity
    fin = open( tmpfile, "r" )
    W = 0.
    nch = 0
    dv = 0.   # velocity width/chan will be computed from v1 and v2
    v1 = -999999.
    v2 = -999999.
    for line in fin :
      a = line.split()
      if len(a) > 2 :
        vlsr = float( a[1] )
        if (vlsr >= vrange[0]) and (vlsr <= vrange[1]) :
          if nch == 0 :
            v1 = vlsr   # v1 is filled in once, is central velocity of first chan
          v2 = vlsr      # v2 is updated each time, will be central velocity of last chan
          nch = nch + 1
          W = W + float( a[3] )
    if (nch > 1) :
      dv = abs(v2 - v1)/(nch - 1)
      print "v1 = %.1f, v2 = %.2f, dv = %.2f km/sec, %d channels in range" % (v1, v2, dv, nch)
    fin.close()
    return [W * dv, nch]
   
     
# set up to plot 5 panels wide by 8 panels high on 8 x 11 sheet, but use only nrows x ncols in upper left corner of sheet;
def plotSnippets2( infile='snip', nrows=3, ncols=2 ) :

  nrowpage = 6
  ncolpage = 4

  Lines = []   # list of Line dictionaries
  fin = open( infile, "rb" )
  while (True) :
    try :
      Line = pickle.load( fin )
      Lines.append( Line )
    except :
      break 
  print "read in %d Line objects" % ( len(Lines) )

  if plotParams["pdf"] :
    pyplot.ioff()
    pp = PdfPages("snippets.pdf")
  ymin = plotParams["ymin"]
  ymax = plotParams["ymax"]
  fig = pyplot.figure( figsize=(8,11) )

# figure out list of positions that will be used
  npList = []
  for np in range(1,nrowpage*ncolpage+1,1) :
    if (np%ncolpage > 0) and (np%ncolpage <= ncols) :
      npList.append( np )
  nn = 0

  for Line in Lines :
    vel = numpy.array( Line["vel"] )
    flux = numpy.array( Line["flux"] )
    if nn > len(npList) :
      if plotParams["pdf"] :
        pyplot.savefig( pp, format='pdf' )
      pyplot.show()
      nn = 0
    p = pyplot.subplot(nrowpage,ncolpage,npList[nn])
    # ymin,ymax = yminmax( numpy.array(Line["flux"]) )
    ymin = -.1 * Line["contLevel"]
    ymax = 2.1*Line["contLevel"]
    plotParams["ymin"] = ymin
    plotParams["ymax"] = ymax
    plotParams["contLevel"] = Line["contLevel"]
    vmin = -60.
    vmax = 70.
    p.axis( [vmin, vmax, ymin, ymax], fontsize=6 )
    p.grid( True, linewidth=.01, color="0.8" )   # bigger color numbers give lighter grays!
    p.plot( vel, flux, color="r", linewidth=1. )
    p.axhline( y=Line["contLevel"], linestyle='--', color='blue')
    p.axhline( y=0, linestyle='-', lw=0.1, color='blue')
    p.tick_params( axis='both', which='major', labelsize=6, length=1 )
    if npList[nn]%ncols == 1 : 
      p.set_ylabel("flux density (Jy)", fontsize=6)
    if npList[nn]-(nrows-1)*ncols > 0  : 
      p.set_xlabel("V$_{LSR}$ (km/sec)", fontsize=6)
  
  # annotate plot, shade the line
    p.text(.04, .9,  "%s" % Line["name"], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
    p.text(.97, .9,  "%.4f GHz" % Line["linefreq"], horizontalalignment='right', transform=p.transAxes, fontsize=6 )
  # comment out following line for recomb line plot
    # p.text(.04, .8,  "T$_L$ = %.0f K" % Line["TL"], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
    #p.text(.97, .8,  "I = %.3f" % Line["intensity"], horizontalalignment='right', transform=p.transAxes)
    yval = yvalue( vel, flux, 5. )
    if yval < Line["contLevel"] : yval = Line["contLevel"]
    print 'yval = %.3f' % yval
    p.annotate( "", xy=(5.,yval), \
             xytext=(5.,0.75*ymax), horizontalalignment="center", rotation="vertical", size=6, \
             arrowprops=dict( arrowstyle="->", linewidth=0.3 ) )
    p.fill_between( vel, Line["contLevel"], flux, color=Line["shadeColor"], alpha=0.5, linewidth=0, \
              where=numpy.abs(vel-5.) < Line["shadeWidth"]/2. )
  # comment in following lines for recomb line plot
    goodRect = Rectangle( (-13.,Line["contLevel"]), 36., 0.2*(ymax-ymin), \
       fill=True, edgecolor='orange',linewidth=0.1, facecolor='yellow', alpha=.5)   # target range
    p.add_patch( goodRect )

  # advance panel number at the very end
    nn = nn + 1
    
  if plotParams["pdf"] :
    pyplot.savefig( pp, format='pdf' )
    pp.close()
  pyplot.show()

# ========================== plot position velocity diagram ========================== #

def extractParameter( instring, param ) :
    '''extract parameter "param" from imlist output, passed to this routine as instring'''
    n = instring.find( param )
    if n > 0 :
      a = instring[n:].split()
      return float( a[2] )
    else :
      return None
  
# ----------------------------------------------------------------------------------------- #
# process pv diagram written out by velplot as a miriad map
# axis1 is VELO; axis2 is position; axis3 is intensity

def readpv( imageFile ) :
    '''read pos-vel map in Miriad format created by velplot, convert to numpy array'''

  # begin by figuring out axes
    p = subprocess.Popen( ( shlex.split('imlist in=%s units=absolute' % imageFile ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    vstart = extractParameter( result, "crval1" )
    pstart = extractParameter( result, "crval2" )
    vstep = extractParameter( result, "cdelt1" )
    pstep = extractParameter( result, "cdelt2" )
    print "vstart = %.2f, vstep = %.2f, pstart = %.3f, pstep = %.3f" % (vstart,vstep,pstart,pstep) 

  # now dump the data and assign to array
    p = subprocess.Popen( ( shlex.split("imtab in=%s log=imtablog format=(3F12.5)" % imageFile ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
    result = p.communicate()[0]
    print result

  # open the logfile, fill out lists, then convert to array
    vel = []
    pos = []
    flx = []
    fin = open("imtablog", "r")
    for line in fin :
      a = line.split()
      vel.append( float(a[0]) + vstart)   # reverse RA to get arcsec on sky
      pos.append( float(a[1]) + pstart)
      flx.append( float(a[2]) )
    fin.close()
    return [ numpy.array(vel), numpy.array(pos), numpy.array(flx) ]

# pv plot
def plotpv( imageFile ) :
    fig = pyplot.figure()
    ax1 = fig.add_axes( [.1,.1,.7,.7])      # main plot
    ax2 = fig.add_axes([.85,.1,.015,.7])    # color wedge
    ax2.tick_params( labelsize=10 )
    ax1.tick_params( labelsize=10 )
    [vel, pos, flx] = readpv( imageFile )
    ax1.axis( [vel[0], vel[-1], pos[0], pos[-1]] )
    nv = 64
    np = 25
    print "nv = %d, np = %d" % (nv,np)
    imgplot = ax1.imshow(numpy.reshape(flx, (np,nv) ), origin='lower', aspect='auto', \
        extent=[vel[0],vel[-1],pos[0],pos[-1]] )
    ax1.grid( True, linewidth=1, color="white")   # color=0.1 is a light gray
    pyplot.colorbar( imgplot, cax=ax2  )
    pyplot.show()

# side-by-side spectra and pv diagrams
# designed to put up to 6 rows x 4 columns on 1 page
def plotSnippets3( infile='snip', nrows=6, ncols=2 ) :

    nrowpage = 6
    ncolpage = 3

    Lines = []   # list of Line dictionaries
    fin = open( infile, "rb" )
    while (True) :
      try :
        Line = pickle.load( fin )
        Lines.append( Line )
      except :
        break 
    print "read in %d Line objects" % ( len(Lines) )
    for Line in Lines:
      print Line["name"],Line["linefreq"]

    pyplot.ioff()
    pp = PdfPages("snippets.pdf")
    fig = pyplot.figure( figsize=(8,11) )

  # make sequential list of positions that will be used
  # panel #1 is upper left corner, panel nrows*ncols is lower right corner

    npList = []
    for np in range(1,nrowpage*ncolpage+1,1) :
      if (np%ncolpage > 0) and (np%ncolpage <= ncols) :
        npList.append( np )
    print "npList = ",npList
    nn = 0  # nn is the index used for npList

  # process only snippets where pvfilename is given
    for Line in Lines :
      vel = numpy.array( Line["vel"] )
      flux = numpy.array( Line["flux"] )

    # if page fills up, save the pdf, display the figure, then reset index nn
      if nn >= len(npList) :
        pyplot.savefig( pp, format='pdf' )
        #pp.close()
        pyplot.show()
        #pp = PdfPages("snippets.pdf")
        fig = pyplot.figure( figsize=(8,11) )
        nn = 0

      #print "nn = %d, npList[nn] = %d" % (nn, npList[nn])
      p = pyplot.subplot(nrowpage,ncolpage,npList[nn])
      ymin,ymax = yminmax( numpy.array(Line["flux"]) )
      #vmin,vmax = yminmax( numpy.array(Line["vel"]) )
         # commented this out because one line is at edge of band, not fully mapped
      vmin = -30.
      vmax = 40.
      ymax = ymax + 0.2*(ymax-ymin)
      plotParams["ymin"] = ymin
      plotParams["ymax"] = ymax
      plotParams["contLevel"] = Line["contLevel"]
      p.axis( [vmin, vmax, ymin, ymax], fontsize=6 )
      p.grid( True, linewidth=.01, color="0.8" )   # bigger color numbers give lighter grays!
      p.plot( vel, flux, color="r", linewidth=1. )
      p.axhline( y=Line["contLevel"], linestyle='--', color='blue')
      p.axhline( y=0, linestyle='-', lw=0.1, color='blue')
      p.tick_params( axis='both', which='major', labelsize=6, length=1 )
      if npList[nn]%ncolpage == 1 : 
        p.set_ylabel("flux density (Jy)", fontsize=6)
      if npList[nn]-(nrows-1)*ncolpage > 0  : 
        p.set_xlabel("V$_{LSR}$ (km/sec)", fontsize=6)

    # annotate plot, shade the line
      p.text(.04, .9,  "%s" % Line["name"], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
      p.text(.97, .9,  "%.4f GHz" % Line["linefreq"], horizontalalignment='right', transform=p.transAxes, fontsize=6 )
      p.fill_between( vel, Line["contLevel"], flux, color=Line["shadeColor"], alpha=0.5, linewidth=0, \
               where=numpy.abs(vel-5.) < Line["shadeWidth"]/2. )
      p.text(.04, .8,  "T$_U$ = %.0f K" % Line["TU"], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
      # p.text(.97, .8,  "I = %.3f" % Line["intensity"], horizontalalignment='right', transform=p.transAxes)

    # this block inserts arrow centered at 5 km/sec; leave it out
      # yval = yvalue( vel, flux, 5. )
      # if yval < Line["contLevel"] : yval = Line["contLevel"]
      # print 'yval = %.3f' % yval
      # p.annotate( "", xy=(5.,yval), \
      #         xytext=(5.,0.75*ymax), horizontalalignment="center", rotation="vertical", size=6, \
      #         arrowprops=dict( arrowstyle="->", linewidth=0.3 ) )

    # advance panel number at the very end
      nn = nn + 1

    # plot pv file if it is available; leave panel blank otherwise
      pvFile = findpvFile( Line["linefreq"] ) 
      print Line["linefreq"], pvFile 
      if pvFile :
        [vel, pos, flx] = readpv( pvFile )
        nv = len( numpy.unique(vel))
        np = len( numpy.unique(pos))
        print "nv = %d, np = %d" % (nv,np)
        p = pyplot.subplot(nrowpage,ncolpage,npList[nn])
        p.tick_params( axis='both', which='major', labelsize=6, length=1 )
        p.axis( [-26., 37., pos[0], pos[-1]] )

      # Deal with situation where there is super-deep absorption next to line
        minflx = numpy.amin( flx )
        maxflx = numpy.amax( flx )
        vminVal = None
        if minflx < -1.*maxflx : vminVal = -1.*maxflx
        imgplot = p.imshow(numpy.reshape(flx, (np,nv) ), origin='lower', aspect='auto', \
            extent=[vel[0],vel[-1],pos[0],pos[-1]], vmin=vminVal )
        p.grid( True, linewidth=0.1, color="white")   # color=0.1 is a light gray
        # pyplot.colorbar( imgplot, cax=ax2  )
      nn = nn + 1
    
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

# this is a special purpose routine to search selected directories for the pv.mp file
#   for spectral line with frequency linefreq
# the prefix to the file can contain any character string, but file must end in
#   "xxx.xx.pv.mp", where xxx.xx is linefreq to 2 decimal places

def findpvFile( linefreq, searchType=".pv.mp" ) :
    dirList = ["/big_scr6/plambeck/350GHz/miriad", "/big_scr6/plambeck/650GHz/miriad"]
      # list of directories that will be searched
    searchstring = "%.2f" % linefreq + searchType
      # filename must include linefreq in GHz, rounded to 2 decimal accuracy
    for directory in dirList :
      for filename in os.listdir( directory ) :
        if searchstring in filename :
          return directory + "/" + filename
    print "couldn't find filename with searchstring %s" % searchstring
    return None
