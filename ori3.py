# ori2.py
# non-OO version of ori.py

import math
import time
import datetime
import dateutil.parser
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
import matplotlib.colors as colors 
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import cm
from scipy import signal
from scipy.stats import norm
import matplotlib.mlab as mlab

ckms = 2.99792e5    # speed of light in km/sec
hplanck = 6.626e-27  # planck's constant in erg-s
kboltzmann = 1.38e-16  # boltzmann's constant in erg/K
ccgs = 2.99792e10   # speed of light in cm/sec

# plotParams dictionary controls plotting details
plotParams = { "npanels" : 2,
               "maxPanelsPerPage" : 2,
               "nlap" : 20,
               "ymin" : .4,       # -.1 for band7, -.2 for band9
               "ymax" : 1.1,      # 1.95 for band7, 4.8 for band9
               #"linelistFile" : "/o/plambeck/OriALMA/Band7B/Spectra/splat_ann.csv", 
               #"linelistFile" : "/o/plambeck/OriALMA/Spectra/splat_ann.csv", 
               "linelistFile" : "/o/plambeck/OriALMA/Salty/LowResData/salts.csv", 
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

#B7spw0 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw0.100plus.a.txt", 
#          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b0"}
#B7spw1 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw1.100plus.a.txt", 
#          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b1"}
#B7spw2 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw2.100plus.a.txt", 
#          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b2"}
#B7spw3 = {'file': "/big_scr6/plambeck/350GHz/miriad/spw3.100plus.a.txt", 
#          'flag': "/big_scr6/plambeck/350GHz/miriad/flags_b3"}

B9spw0 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw0.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a0"}
B9spw1 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw1.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a1"}
B9spw2 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw2.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a2"}
B9spw3 = {'file': "/big_scr6/plambeck/650GHz/miriad/spw3.100plus.a.txt", 
          'flag': "/big_scr6/plambeck/650GHz/miriad/flags_a3"}

BN7spw0 = {'file': "/big_scr6/plambeck/350GHzBN/miriad/spw0.ch100.txt", 
           'flag': "/o/plambeck/OriALMA/Band7/Spectra/BNflags_b0"}
BN7spw1 = {'file': "/big_scr6/plambeck/350GHzBN/miriad/spw1.ch100.txt", 
           'flag': "/o/plambeck/OriALMA/Band7/Spectra/BNflags_b1"}
BN7spw2 = {'file': "/big_scr6/plambeck/350GHzBN/miriad/spw2.ch100.txt", 
           'flag': "/o/plambeck/OriALMA/Band7/Spectra/BNflags_b2"}
BN7spw3 = {'file': "/big_scr6/plambeck/350GHzBN/miriad/spw3.ch100.txt", 
           'flag': "/o/plambeck/OriALMA/Band7/Spectra/BNflags_b3"}

B8spw0 = {'file': "/big_scr6/plambeck/460GHz/miriad/spw0.ch250.txt", 
           'flag': "/o/plambeck/OriALMA/Band8/Spectra/flags_a0"}
B8spw1 = {'file': "/big_scr6/plambeck/460GHz/miriad/spw1.ch250.txt", 
           'flag': "/o/plambeck/OriALMA/Band8/Spectra/flags_a1"}
B8spw2 = {'file': "/big_scr6/plambeck/460GHz/miriad/spw2.ch250.txt", 
           'flag': "/o/plambeck/OriALMA/Band8/Spectra/flags_a2"}
B8spw3 = {'file': "/big_scr6/plambeck/460GHz/miriad/spw3.ch250.txt", 
           'flag': "/o/plambeck/OriALMA/Band8/Spectra/flags_a3"}

B6spw0 = {'file': "/big_scr6/plambeck/220GHz/miriad/spw0.ch100.txt", 
           'flag': "/big_scr6/plambeck/220GHz/miriad/flags_a0"}
B6spw1 = {'file': "/big_scr6/plambeck/220GHz/miriad/spw1.ch100.txt", 
           'flag': "/big_scr6/plambeck/220GHz/miriad/flags_a1"}
B6spw2 = {'file': "/big_scr6/plambeck/220GHz/miriad/spw2.ch100.txt", 
           'flag': "/big_scr6/plambeck/220GHz/miriad/flags_a2"}
B6spw3 = {'file': "/big_scr6/plambeck/220GHz/miriad/spw3.ch100.txt", 
           'flag': "/big_scr6/plambeck/220GHz/miriad/flags_a3"}

B7Bspw0 = {'file': "/alma_scr/plambeck/340GHz/miriad/spw0.ch100.txt", 
           'flag': "/alma_scr/plambeck/340GHz/miriad/flags_a0"}
B7Bspw1 = {'file': "/alma_scr/plambeck/340GHz/miriad/spw1.ch100.txt", 
           'flag': "/alma_scr/plambeck/340GHz/miriad/flags_a1"}
B7Bspw2 = {'file': "/alma_scr/plambeck/340GHz/miriad/spw2.ch100.txt", 
           'flag': "/alma_scr/plambeck/340GHz/miriad/flags_a2"}
B7Bspw3 = {'file': "/alma_scr/plambeck/340GHz/miriad/spw3.ch100.txt", 
           'flag': "/alma_scr/plambeck/340GHz/miriad/flags_a3"}

# put these in on 17oct2018 to run from harpo on /o/plambeck/OriALMA/Salty/LowResData
B7spw0 = {'file': "spw0.100plus.a.txt", 
          'flag': "flags_b0"}
B7spw1 = {'file': "spw1.100plus.a.txt", 
          'flag': "flags_b1"}
B7spw2 = {'file': "spw2.100plus.a.txt", 
          'flag': "flags_b2"}
B7spw3 = {'file': "spw3.100plus.a.txt", 
          'flag': "flags_b3"}

specList = [B7spw0,B7spw1,B7spw2,B7spw3,B9spw0,B9spw1,B9spw2,B9spw3]


# run getspec, dump spectrum to a text file; adjust freq scale to be correct for source at LSR velocity vsource
def dumpspec( infile, outfile, region='arcsec,box(.1)', vsource=0., hann=1 ):
    [chan, freq, flux ] = getspec( infile, region=region, vsource=vsource, hann=hann )
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
    vs0 = 0.    # default value
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

  # special case: renumber channels to handle B3 BN spectra made in 4 sections, each labeled chan 1-480
    if len(chan) == 1920 and chan[0] != 1920 and chan[-1] != 1920 :
      if chan[0] == 1 :
        chan = range(1,1921,1)
      else :
        chan = range(1920,0,-1)
    return numpy.array(chan), numpy.array(freq), numpy.array(flux)

# ----------------------------------------------------------------------------------------------
# getspec returns spectrum of selected region from spectral line cube; note that region
#  command could also select channel range, if desired
# Note: - axis3 of spectral line data cube is LSR velocity relative to restfreq (at LSR=0)
#       - from this, one can compute the restfreq of that channel in the LSR=0 frame (the lab frame) 
#       - getspec divides by doppler factor to compute restfreq in the source frame (lsr=vsource) 

def getspec( infile, region='relpix,box(-2,-2,0,0)', vsource=0., hann=1, tmpfile="junk" ):
    '''dump out spectrum of selected region with imspec, return [chan, freqLSR, flux] arrays'''

  # step 1: use imlist to retrieve velocity and freq information from the header
    restfreq,v1,dv = getVelInfo( IQUVmapList[0] )
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

def getuvspec( infile, vsource=5., hann=3, tmpfile="junk", select='uvrange(0,200)', offset="0.0,0.0" ):
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

  # CAUTION: this used to use option=ampscalar!
    p= subprocess.Popen( ( shlex.split("uvspec vis=%s select=%s axis=freq,amp stokes=i \
      options=avall,nobase interval=100000 hann=%d log=%s offset=%s" % \
      (infile,select,hann,tmpfile,offset) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)

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

def processLineList( lineListFile, vsource ) :
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
          name.append( a[nameCol] )
          vdop.append( vlsr )
          dopfac = 1. - (vlsr - vsource)/ckms
            # ... doppler correct lab freq if line is emitted at velocity different than vsource

        # sometimes freq shows up as calculated freq, other times as measured freq
          if len(a[freqCol]) > 0 :
            freq.append( dopfac * float(a[freqCol].split(",")[0])/1.e9 )
            print "... ", a[nameCol], dopfac, a[freqCol], freq[-1]
            #freq.append( dopfac * float(a[freqCol]) )
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
# parse splatalogue output, remove unobserved frequencies, duplicates u

#def pruneLineList( infile="/o/plambeck/Downloads/splatalogue.csv", outfile="/o/plambeck/OriALMA/Spectra/extra.csv" ) :
#def pruneLineList( infile="/o/plambeck/OriALMA/Band6/Spectra/splatalogue.csv", outfile="/o/plambeck/OriALMA/Band6/Spectra/splat_ann.csv" ) :
def pruneLineList( infile="/o/plambeck/Downloads/splatalogue.csv", outfile="/o/plambeck/OriALMA/Band7B/Spectra/splat_ann.csv" ) :
    #freqRanges = [ [340.545, 344.310], [352.665, 356.430], [649.282, 651.177], [661.303, 665.067], [665.922, 667.817] ]
    #freqRanges = [ [462.5, 463.5], [463.6, 465.6], [473.9, 474.9], [475.7, 477.7] ]
    #freqRanges = [ [229.17, 231.05], [231.83, 233.71], [214.28, 216.16], [216.98, 218.86] ]
    freqRanges = [ [344.05, 345.93], [346.05,347.93], [333.94,335.82], [332.05, 333.93] ]
    copy = []
    name = []
    QN = []
    freq = []
    cat = []
    EL = []

  # pass 1: read file, create lists of species, freq, QN, linelist for each input line 
    fin = open( infile, "r" )
    for line in fin :
      copy.append( 1 )     # default is to copy this line
      name.append( "" )      
      freq.append( 0. )
      QN.append( "" )
      cat.append( "" )
      EL.append( 0. )
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
            if "E_L (K)" in a[n] :
              ELCol = n
              print "ELCol = %d" % ELCol
        else : 
          name[-1] = a[nameCol]
          freq[-1] = float(a[freqCol])
          QN[-1] = a[QNcol]
          cat[-1] = a[catCol]
          EL[-1] = float(a[ELCol])

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
            # print freq[m],cat[m]
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
        if EL[n] > 500. :
        #fout.write("%s:None:None:\n" % line.rstrip('\n'))
          fout.write("%s:yellow:36.:\n" % line.rstrip('\n'))
        else :
          fout.write("%s:green:36.:\n" % line.rstrip('\n'))
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
def ann( p, freq, flux, fmin, fmax, plotParams, vsource ) :
    '''annotates plot with spectral line labels'''
    print "annotating freq range %.3f to %.3f GHz" % (fmin,fmax)
    if plotParams["linelistFile"] :
      name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( plotParams["linelistFile"], vsource )
      for n in range(0,len(name)) :
        if (linefreq[n] >= fmin) and (linefreq[n] <= fmax) :
          print "fmin = %.4f, fmax = %.4f, linefreq = %.4f" % (fmin,fmax,linefreq[n])
          yval = yvalue( freq, flux, linefreq[n] )
          if yval < plotParams["contLevel"] : yval = plotParams["contLevel"]
          if yval > 0.8*plotParams["ymax"] : yval = plotParams["contLevel"]

          #if TU[n] > 0. :
          #  p.annotate( "%s  %.0fK" % (name[n],TU[n]), xy=(linefreq[n],yval), \
          #   xytext=(linefreq[n],0.92*plotParams["ymax"]), horizontalalignment="center", rotation="vertical", size=7, \
          #   arrowprops=dict( arrowstyle="-", linewidth=0.2 ) )    # use arrowstyle="->" to get arrow
          #else :
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

def findRMS( chan, freq, flux, flagtable, nchans=1920 ) :

    flaggedBad = numpy.zeros( len(chan), dtype=int )
    if flagtable :
      try :
        fin = open( flagtable, "r" )
      except :
        print "couldn't find flagtable %s - skipping rms and contLevel calculation" % flagtable
        return [0.,0.,flaggedBad] 

    # note that flagtable is just a list of channels "nchans" in the original uv data;
    # there may be a different number of channels "nch" in the map if it was made
    #    with line=chan,nch,1,nstep,nstep
    # figure out nstep based on nchans (input parameter) and nch (length of array) 
    # for B7 and B9 data, nchans=3840 (the default), but for B8 Hirota data, nchans=1920
    # figure out nscale from number of channels in spectrum; then, to prevent mistakes,
    #   confirm that frequency increment makes sense
    # note that there is no protection against mixing up flagtables
      
      
      nscale = nchans/len(chan)
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
        ann( p, freq, fluxList[0], freq[nch1], freq[nch2], plotParams, vsource )

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

# CAUTION: edit nchans in findRMS

#def pubFig( specList=[B9spw0,B9spw1,B9spw2,B9spw3], plotParams=plotParams, vsource=5. ) :
#def pubFig( specList=[B7spw0,B7spw1,B7spw2,B7spw3], plotParams=plotParams, vsource=5. ) :
#def pubFig( specList=[BN7spw0,BN7spw1,BN7spw2,BN7spw3], plotParams=plotParams, vsource=9. ) :
#def pubFig( specList=[BN7spw0,BN7spw1], plotParams=plotParams, vsource=9. ) :
#def pubFig( specList=[B8spw0,B8spw1,B8spw2,B8spw3], plotParams=plotParams, vsource=5. ) :
#def pubFig( specList=[B6spw0,B6spw1,B6spw2,B6spw3], plotParams=plotParams, vsource=5. ) :
#def pubFig( specList=[B7Bspw0,B7Bspw1,B7Bspw2,B7Bspw3], plotParams=plotParams, vsource=5. ) :

def pubFig( specList, plotParams, vsource=5. ) :
    pyplot.ioff()
    pp = PdfPages("spectrum.pdf")
    ymin = plotParams["ymin"]
    ymax = plotParams["ymax"]
    fig = pyplot.figure( figsize=(8,11) )

    # fig = pyplot.figure( figsize=(11,8) )

    npanel = 1
    for spectrum in specList: 
      [chan, freq, flux ] = readspec( spectrum["file"], vsource=vsource )
      p = pyplot.subplot( plotParams["npanels"], 1, npanel)
      p.grid( True, linewidth=0.05, color="0.01" )   # color=0.1 is a light gray

      fmin = freq[0]
      fmax = freq[-1]
      delta = fmax - fmin
      fmin = fmin - .04*delta
      fmax = fmax + .04*delta

      chmin = chan[0]
      chmax = chan[-1]
      delta = chmax - chmin    # note that delta can be negative!
      chmin = chmin - .04*delta
      chmax = chmax + .04*delta

      p.axis( [fmin, fmax, ymin, ymax] )
      p3 = p.twiny()   
      p3.axis( [chmin, chmax, ymin, ymax] )
      #p3.xaxis.set_minor_locator( x_locator )
      p3.tick_params( axis='x', which='major', labelsize=8 )

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
        ann( p, freq, flux, freq[0], freq[-1], plotParams, vsource )
      if npanel == plotParams["npanels"] :
         p.set_xlabel("freq (GHz)", fontsize=8)
      p.set_ylabel("flux density (Jy)", fontsize=8)
      if npanel >= plotParams["maxPanelsPerPage"] :
        pyplot.savefig( pp, format='pdf' )
        pyplot.show()
        npanel = 1
      else :
        npanel = npanel + 1

    if npanel > 1 :
      pyplot.savefig( pp, format='pdf' )
      pyplot.show()
    pp.close()

# returns line intensities in LINEAR units for a particular temperature T
def calcIntensity( infile, T ) :
  T0 = 300.
  hPlanck = 6.62607e-27 
  kBoltzmann = 1.38065e-16

  name,freq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( infile, vsource=0. )
  Iba = numpy.zeros( len(name) )
  for n in range(0, len(name) ) :
    TU = TL[n] + (hPlanck * freq[n]*1.e9 / kBoltzmann)   # upper state energy in K
    Iba[n] = pow(10.,intensity[n]) * \
      (math.exp(-TU/T) - math.exp(-TL[n]/T)) / (math.exp(-TU/T0) - math.exp(-TL[n]/T0))
  return Iba

def plotIntensity( infile ) :
  name,freq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( infile, vsource=0. )
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
    vrange=90., vsource=5.0 ) :

  # read the lineList file
    name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( lineListFile, vsource )

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
                print "red wing :  %.3f  %.3f" % (vel, flux[i]-Line["contLevel"] )
              if (vel > -13.) and (vel <= 5.) :
                Line["intfluxBlue"] = Line["intfluxBlue"] + (flux[i] - contLevel)
                nblue = nblue + 1
                print "blue wing:  %.3f  %.3f" % (vel, flux[i]-Line["contLevel"] )

            else :     # past end of snippet
              if FillingIn :
                Line["rmsRed"] = rms * math.sqrt(nred/2.)     # assume that spectra were Hanning smoothed! 
                Line["rmsBlue"] = rms * math.sqrt(nblue/2.)   # assume that spectra were Hanning smoothed! 
                print "intRed,Blue = ",Line["intfluxRed"],Line["intfluxBlue"]
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
csvFileList = [ '90GHz.csv', '229GHz.csv', '350GHz.csv', '460GHz.csv', '660GHz.csv' ]
# csvFileList = [ '6.1GHz.csv', '229GHz_500.csv', '350GHz_500.csv', '660GHz_500.csv' ]
srcNameList = [ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l' ]
def plotFluxes( ) :
  nx = 3
  ny = 2
  pyplot.ioff()
  pp = PdfPages("pointSrcFluxes.pdf")
  nplot=0
  for srcName in srcNameList :
    nplot = nplot+1
    p = pyplot.subplot(ny,nx,nplot)
    p.axis( [70., 1000., 1., 10000.] )
    p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
    frqGHz = []
    melfrqGHz = []
    SmJy = []
    uncSmJy = []
    melJy = []
    print " "
    for csvFile in csvFileList :
      fin = open( csvFile, 'r' )
      for line in fin :
        if len(line) > 0 :
          a = line.split(',')
          if (len(a[0]) > 0 ) and (a[0] == srcName) :
          # values from me
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
          # values from Mel
            try: 
              mJy = 1000. * float(a[22])
              if mJy > 0. :
                n = csvFile.find("GHz") 
                melfrqGHz.append( float( csvFile[0:n] ) )
                melJy.append( mJy )
              else :
                pass     # don't include if flux is zero
            except :             # handles case where flux = '-'
              pass
              
    print "source %s" % srcName
    print "frqGHz, SmJy, uncSmJy: ",frqGHz,SmJy,uncSmJy
    print "melGHz, melJy: ", melfrqGHz, melJy

    if len(frqGHz) == len(SmJy) :
      Sf2 = []
      for frq,S in zip(frqGHz,SmJy) : 
        if frq > 200. and frq < 240. :
          fref = frq
          Sref = S
      for frq,S in zip(frqGHz,SmJy) : 
          Sf2.append( pow(frq/fref, 2.) * Sref )
      p.loglog( frqGHz, SmJy )
      p.loglog( melfrqGHz, melJy, color='b', marker='o') 
      p.loglog( frqGHz, Sf2, linestyle='dashed')
      p.errorbar( frqGHz, SmJy, yerr=uncSmJy)
      p.text( .1, .9, "%s" % srcName, transform=p.transAxes )
      p.tick_params( axis='both', which='major', labelsize=8 )
    else :
      print srcName, frqGHz, SmJy
    if nplot >= (nx*ny) :
      pyplot.tight_layout( )
      pyplot.savefig( pp, format='pdf' )
      pyplot.show()
      nplot = 0
      
  if nplot > 0 :
    pyplot.tight_layout( )
    pyplot.savefig( pp, format='pdf' )
    pyplot.show()
  pp.close()

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
BNbox = 'arcsec,box(.14,-.05,-.16,.25)'

def makespec220() :
  dumpspec( "spw0.ch100.cm", "spw0.ch100.txt", region="arcsec,box(.2)", vsource=0., hann=3 )
  dumpspec( "spw1.ch100.cm", "spw1.ch100.txt", region="arcsec,box(.2)", vsource=0., hann=3 )
  dumpspec( "spw2.ch100.cm", "spw2.ch100.txt", region="arcsec,box(.2)", vsource=0., hann=3 )
  dumpspec( "spw3.ch100.cm", "spw3.ch100.txt", region="arcsec,box(.2)", vsource=0., hann=3 )

def makespec460() :
  dumpspec( "spw0.ch250.cm", "spw0.ch250.txt", region="arcsec,box(.2)", vsource=0., hann=3 )
  dumpspec( "spw1.ch250.cm", "spw1.ch250.txt", region="arcsec,box(.2)", vsource=0., hann=3 )
  dumpspec( "spw2.ch250.cm", "spw2.ch250.txt", region="arcsec,box(.2)", vsource=0., hann=3 )
  dumpspec( "spw3.ch250.cm", "spw3.ch250.txt", region="arcsec,box(.2)", vsource=0., hann=3 )

def makespec() :
  dumpspec( "spw0.100plus.cm", "spw0.100plus.c.txt", region=regionc, vsource=0., hann=3 )
  dumpspec( "spw1.100plus.cm", "spw1.100plus.c.txt", region=regionc, vsource=0., hann=3 )
  dumpspec( "spw2.100plus.cm", "spw2.100plus.c.txt", region=regionc, vsource=0., hann=3 )
  dumpspec( "spw3.100plus.cm", "spw3.100plus.c.txt", region=regionc, vsource=0., hann=3 )

# BNmakespec created 1 feb 2017
def BNmakespec() :
  dumpspec( "spw0.ch500.cm", "spw0.ch500.txt", region=BNbox, vsource=0., hann=3 )
  dumpspec( "spw1.ch500.cm", "spw1.ch500.txt", region=BNbox, vsource=0., hann=3 )
  #dumpspec( "spw0.ch100.cm", "spw0.ch100.txt", region=BNbox, vsource=0., hann=3 )
  #dumpspec( "spw1.ch100.cm", "spw1.ch100.txt", region=BNbox, vsource=0., hann=3 )
  #dumpspec( "spw2.ch100.cm", "spw2.ch100.txt", region=BNbox, vsource=0., hann=3 )
  #dumpspec( "spw3.ch100.cm", "spw3.ch100.txt", region=BNbox, vsource=0., hann=3 )

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
def rotdiag2( inlistFile="/o/plambeck/OriALMA/Spectra/snip_rotdiag.csv", outfile="rotdiagpiece.wip", vrange=[5.,20.], tablefile="rotdiagpiece.tab", \
       summaryfile="rotdiagsummary.txt" ) :

    kB = 1.38e-16    # Boltzman's constant in cgs units

  # read list of lines that must be processed
    name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( inlistFile )

  # open output files file, write out header info
    fout = open( outfile, "w" )        # wip commands
    fout.write("# created using ori3.rotdiag2(); using velocity range %.1f to %.1f\n" % (vrange[0],vrange[1]))
    fout2 = open( tablefile, "w" )     # latex table
    fout2.write("# created using ori3.rotdiag2(); using velocity range %.1f to %.1f\n" % (vrange[0],vrange[1]))
    fout3 = open( summaryfile, "w")    # text file with results
    fout3.write("# created using ori3.rotidag2(); using velocity range %.1f to %.1f\n" % (vrange[0],vrange[1]))
    fout3.write("# freq  Sdv   unc   logS   loguncS   Tupper\n")

  # process lines one at a time; find "mc" map, compute integrated intensity
    for n in range(0,len(linefreq)) :
      mcFile = findpvFile( linefreq[n], searchType=".mc" ) 
      if mcFile :
        Sdv,nch,vel,TB = getW( mcFile, vrange )   # get integrated emission; units must be K-km/sec

      # uncertainty estimates are based on BLANK_354.25.mc (~1 K rms) and BLANK_649.80.mc (~ 3K rms);
      #   assume that uncertainty in integrated Tb is sqrt(nch) * unc/chan
        rms = 1.5
        if linefreq[n] > 500. : rms = 5.
        unc = math.sqrt(nch) * rms
        print "W = %.2f unc = %.2f" % (Sdv,unc)
        const = 3.* kB * 1.e5/(8. * pow(math.pi,3) * linefreq[n] * 1.e9 * Smusq[n] * pow(1.e-18, 2) )
          # formula A4 in Cummins et al 1986; 1.e5 is to convert km/sec to cm/sec; 1.e-36 converts debye^2 to (esu-cm)^2
        fout3.write(" %11.5f  %10.2f %7.2f  %10.4f  %10.4f   %8.1f\n" % \
         ( linefreq[n], Sdv, unc, evalLog(const*Sdv), evalLog(const*Sdv) - evalLog(const*(Sdv-unc)), TU[n] ))
        print " %11.5f  %10.2f %7.2f  %10.4f  %10.4f  %10.4f  %8.1f" % \
         ( linefreq[n], Sdv, unc, evalLog(const*Sdv), evalLog(const*(Sdv-unc)), evalLog(const*(Sdv+unc)), TU[n] )

      # fout2 output is used for table in ms.tex
        fout2.write("%4.0f &  %11.3f & %6.2f & %5.0f (%.0f) \\\\ \n" % \
         ( TU[n], linefreq[n], Smusq[n], Sdv, unc) )
        fout.write("# %11.5f  %10.2f %7.2f  %10.4f  %8.1f  %8.3f  %s  %s\n" % \
         ( linefreq[n], Sdv, unc, evalLog(const*Sdv),TU[n],Smusq[n],name[n],QN[n]) )
        fout.write("move %.2f %.4e\n" % (TU[n], evalLog(const*Sdv) ) )
        fout.write("dot\n")
        fout.write("move %.2f %.4e\n" % (TU[n], evalLog(const*(Sdv-unc)) ) )
        fout.write("draw %.2f %.4e\n" % (TU[n], evalLog(const*(Sdv+unc)) ) )
      else :
        print "skipping line at freq %.3f - couldn't find file %s" % (linefreq[n],mcFile)
    fout.close()
    fout2.close()
    fout3.close()

# helper routine that opens "mc" map file (line minus continuum), uses imspec to compute mean Tb in region
# revised 18 apr to return velocity, TB arrays as well
def getW( mcFile, vrange, region="arcsec,box(0.1)", tmpfile="junk" ) :

    vel = []
    TB = []
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
        vel.append(vlsr)
        TB.append( float(a[3] ))
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
    return [W * dv, nch, numpy.array(vel), numpy.array(TB) ]
   
     
# set up to plot ncolpage panels wide by nrowpage panels high on 8 x 11 sheet, 
#   but use only nrows x ncols in upper left corner of sheet;
def plotSnippets2( infile='snip', nrows=4, ncols=2, vmin=-40., vmax=50. ) :

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

  pyplot.ioff()
  pp = PdfPages("snippets.pdf")
  ymin = plotParams["ymin"]
  ymax = plotParams["ymax"]
  fig = pyplot.figure( figsize=(8,11) )

# figure out list of positions that will be used
  npList = []
  for np in range(1,nrows*ncolpage+1,1) :
    if (np%ncolpage > 0) and (np%ncolpage <= ncols)  :
      npList.append( np )
  print npList
  nn = 0

  for Line in Lines :
    vel = numpy.array( Line["vel"] )
    flux = numpy.array( Line["flux"] )
    p = pyplot.subplot(nrowpage,ncolpage,npList[nn])
    # ymin,ymax = yminmax( numpy.array(Line["flux"]) )
    ymin = -.1 * Line["contLevel"]
    ymax = 2.1*Line["contLevel"]
    plotParams["ymin"] = ymin
    plotParams["ymax"] = ymax
    plotParams["contLevel"] = Line["contLevel"]
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
    p.text(.04, .8,  "E$_U$ = %.0f K" % Line["TU"], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
    #p.text(.97, .8,  "I = %.3f" % Line["intensity"], horizontalalignment='right', transform=p.transAxes)
    yval = yvalue( vel, flux, 5. )
    if yval < Line["contLevel"] : yval = Line["contLevel"]
    print 'yval = %.3f' % yval
    p.annotate( "", xy=(5.,yval), \
             xytext=(5.,0.75*ymax), horizontalalignment="center", rotation="vertical", size=6, \
             arrowprops=dict( arrowstyle="->", linewidth=0.3 ) )
    p.fill_between( vel, Line["contLevel"], flux, color=Line["shadeColor"], alpha=0.5, linewidth=0, \
              where=numpy.abs(vel-5.) < Line["shadeWidth"]/2. )

  # comment in following lines for recomb line plot; plot shaded box over presumed location of line
    # goodRect = Rectangle( (-13.,Line["contLevel"]), 36., 0.2*(ymax-ymin), \
    #    fill=True, edgecolor='orange',linewidth=0.1, facecolor='yellow', alpha=.5)   # target range
    # p.add_patch( goodRect )

  # advance panel number at the very end
    nn = nn + 1
# temporary kludge to avoid printing empty page
    #if nn >= len(npList) :
    #  pyplot.savefig( pp, format='pdf' )
    #  pyplot.show()
    #  nn = 0
    print "nn = %d" % nn
    
  pyplot.savefig( pp, format='pdf' )
  pp.close()
  pyplot.show()

# copy of plotSnippets2 specially tailored for recomb lines and overlay of SO2
def plotRecomb( infile='snip_recomb', nrows=2, ncols=1, vmin=-45., vmax=55. ) :

  matplotlib.rcParams['axes.linewidth'] = 0.7
    # GLOBALLY set width of plot borders - better exit and restart after using this routine!

  nrowpage = 5
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

  pyplot.ioff()
  pp = PdfPages("recomb3.pdf")
  ymin = plotParams["ymin"]
  ymax = plotParams["ymax"]
  fig = pyplot.figure( figsize=(8,11) )

# figure out list of positions that will be used
  npList = [1,4,4,4]
  print npList
  nn = 0

  for Line in Lines :
    vel = numpy.array( Line["vel"] )
    flux = numpy.array( Line["flux"] )
    p = pyplot.subplot(nrowpage,ncolpage,npList[nn])
    ymin = -.1 * Line["contLevel"]
    ymax = 2.1*Line["contLevel"]
    plotParams["ymin"] = ymin
    plotParams["ymax"] = ymax
    plotParams["contLevel"] = Line["contLevel"]
    p.axis( [vmin, vmax, ymin, ymax], fontsize=6, lwidth=0.5 )
    p.grid( True, linewidth=.01, color="0.8" )   # bigger color numbers give lighter grays!
    if nn == 2 :
      nstart = 0
      nstop = 0
      for v in vel :
        if v > 30. : nstart = nstart + 1
        if v > -20. : nstop = nstop + 1
      print vel[nstart:nstop]
      p.plot( vel[nstart:nstop], flux[nstart:nstop], color="g", linewidth=0.5 )
      p.text(.04, .84,  "%s" % Line["name"], horizontalalignment='left', transform=p.transAxes, fontsize=6, color='g' )
      p.text(.97, .84,  "%.4f GHz" % Line["linefreq"], horizontalalignment='right', transform=p.transAxes, fontsize=6, color='g' )
    elif nn == 3 :
      nstart = 0
      nstop = 0
      for v in vel :
        if v > 30. : nstart = nstart + 1
        if v > -20. : nstop = nstop + 1
      print vel[nstart:nstop]
      p.plot( vel[nstart:nstop], flux[nstart:nstop], color="b", linewidth=0.5  )
      p.text(.04, .78,  "%s" % Line["name"], horizontalalignment='left', transform=p.transAxes, fontsize=6, color='b' )
      p.text(.97, .78,  "%.4f GHz" % Line["linefreq"], horizontalalignment='right', transform=p.transAxes, fontsize=6, color='b' )
    else :
      p.plot( vel, flux, color="red", linewidth=1 )
      p.plot( [7.,21.], [1.05*Line["contLevel"],1.05*Line["contLevel"]], color="black", linewidth=0.6 )
      p.text(14., 1.07*Line["contLevel"], "limit", horizontalalignment='center', fontsize=6, color='black' )
      p.axhline( y=Line["contLevel"], linestyle='--', color='blue')
      p.axhline( y=0, linestyle='-', lw=0.1, color='blue')
      p.tick_params( axis='both', which='major', labelsize=6, length=1 )
      p.text(.04, .9,  "%s" % Line["name"], horizontalalignment='left', transform=p.transAxes, fontsize=6, color='red' )
      p.text(.97, .9,  "%.4f GHz" % Line["linefreq"], horizontalalignment='right', transform=p.transAxes, fontsize=6, color='red' )
      p.set_ylabel("flux density (Jy)", fontsize=6)
      yval = Line["contLevel"]
      print 'yval = %.3f' % yval
      p.annotate( "", xy=(5.,yval), \
         xytext=(5.,0.75*ymax), horizontalalignment="center", rotation="vertical", size=6, \
         arrowprops=dict( arrowstyle="->", linewidth=0.3 ) )
    if npList[nn]-(nrows-1)*ncols > 0  : 
      p.set_xlabel("V$_{LSR}$ (km/sec)", fontsize=6)
  
  # comment in following lines for recomb line plot; plot shaded box over presumed location of line
    if (nn == 0) or (nn == 3) :
      goodRect2 = Rectangle( (5.,Line["contLevel"]), 18., 0.2*(ymax-ymin), \
         fill=True, edgecolor='none', facecolor='yellow', alpha=.5)   # target range
      goodRect1 = Rectangle( (-13.,Line["contLevel"]), 36., 0.2*(ymax-ymin), \
         fill=False, edgecolor='black',linewidth=0.5, facecolor='yellow', alpha=.5)   # target range
      p.add_patch( goodRect2 )
      p.add_patch( goodRect1 )

  # advance panel number at the very end
    nn = nn + 1
    print "nn = %d" % nn
    
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
# 4/20/16 added vmin,vmax,pmin,pmax as a quick way of selecting piece of image that I want;
#   hopefully defaults are so large that they won't trip me up in cases where I don't care

def readpv( imageFile, vmin=-1.e8, vmax=1.e8, pmin=-1.e8, pmax=1.e8 ) :
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
      v = float(a[0]) + vstart
      p = float(a[1]) + pstart
      if (v >= vmin) and (v <= vmax) and (p >= pmin) and (p <= pmax) :
        vel.append( v )   
        pos.append( p )
        flx.append( float(a[2]) )
    fin.close()
    return [ numpy.array(vel), numpy.array(pos), numpy.array(flx) ]

# pv plot
def plotpv( imageFile, pmin=-1.1, pmax=1.5 ) :
    fig = pyplot.figure()
    ax1 = fig.add_axes( [.1,.1,.7,.7])      # main plot
    ax2 = fig.add_axes([.85,.1,.015,.7])    # color wedge
    ax2.tick_params( labelsize=10 )
    ax1.tick_params( labelsize=10 )
    [vel, pos, flx] = readpv( imageFile, pmin=pmin, pmax=pmax )
    ax1.axis( [vel[0], vel[-1], pos[0], pos[-1]] )
    nv = len( numpy.unique( vel ) )
    np = len( numpy.unique( pos ) )
    print "nv = %d, np = %d" % (nv,np)
    imgplot = ax1.imshow(numpy.reshape(flx, (np,nv) ), origin='lower', aspect='auto', \
        extent=[vel[0],vel[-1],pos[-1],pos[0]], vmin=-.02, vmax=.1 )
    ax1.grid( True, linewidth=1, color="white")   # color=0.1 is a light gray
    pyplot.colorbar( imgplot, cax=ax2  )
    pyplot.show()

# side-by-side spectra and pv diagrams; used to make snip_pv.pdf
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
      vmin = -29.
      vmax = 39.
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
      if (Line["TU"] > 0.) :
        p.text(.04, .8,  "E = %.0f K" % Line["TU"], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
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
        [vel, pos, flx] = readpv( pvFile, -19., 29., -.3, .3 )
            # note: extra parameters limit velocity and position range that is read in
        nv = len( numpy.unique(vel))
        np = len( numpy.unique(pos))
        print "nv = %d, np = %d" % (nv,np)
        p = pyplot.subplot(nrowpage,ncolpage,npList[nn])
        p.tick_params( axis='both', which='major', labelsize=6, length=1 )

      # for this figure, flip orientation of yaxis so SE is at bottom
        #p.axis( [-26., 37., 0.3, -0.3] )
        p.axis( [vel[0],vel[-1],pos[-1],pos[0]] )
        p.set_ylabel("offset (arcsec)", fontsize=6)
        if npList[nn]-(nrows-1)*ncolpage > 0  : 
          p.set_xlabel("V$_{LSR}$ (km/sec)", fontsize=6)

      # Deal with situation where there is super-deep absorption next to line
      # I have now replaced that - just plot for 0 to max
        # minflx = numpy.amin( flx )
        maxflx = numpy.amax( flx )
        vminVal = None
        # if minflx < -1.*maxflx : vminVal = -1.*maxflx
        imgplot = p.imshow(numpy.reshape(flx, (np,nv) ), clim=(0.,maxflx), origin='lower', aspect='auto', \
            extent=[vel[0],vel[-1],pos[0],pos[-1]], vmin=vminVal )
        p.grid( True, linewidth=0.1, color="white")   # color=0.1 is a light gray
        # pyplot.colorbar( imgplot, cax=ax2  )
      nn = nn + 1
   
  # allow a little more space between columns for ylabel
    pyplot.subplots_adjust(wspace=0.3) 
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

# read splat_ann; sort lines as desired; then write out in latex form
def makeLineTable( linelistFile="splat_ann_table.csv", outfile="linetable.tex" ) :
    name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( plotParams["linelistFile"] )
    index = numpy.argsort( linefreq )
    fout = open( outfile, "w" )
    for i in range(0, len(index) ) :
      n = index[i]
      print "%10.3f & %s %s & %5.0f \\" % (linefreq[n],name[n],QN[n],TU[n]) 
      fout.write("%10.4f & %s & %5.0f \\\\\n" % (linefreq[n],name[n],TU[n]) )
    fout.close()
    
# plot SO2 lines used in rotdiag
# set up to plot ncolpage panels wide by nrowpage panels high on 8 x 11 sheet, 
#   but use only nrows x ncols in upper left corner of sheet;
def plotSnippets4( infile='/o/plambeck/OriALMA/Spectra/snip_SO2rotdiag.csv', nrows=8, ncols=2, vmin=-30., vmax=40. ) :

    nrowpage = 8
    ncolpage = 4
    pyplot.ioff()
    pp = PdfPages("snippets.pdf")
    fig = pyplot.figure( figsize=(8,11) )

  # figure out list of positions that will be used
    npList = []
    for np in range(1,nrowpage*ncolpage+1,1) :
      if (np%ncolpage > 0) and (np%ncolpage <= ncols) :
        npList.append( np )
    nn = 0

  # read list of lines that must be processed
    name,linefreq,QN,TL,TU,Smusq,intensity,vdop,shadeColor,shadeWidth = processLineList( infile )

  # process lines one at a time; find "mc" map, compute integrated intensity
    for n in range(0,len(linefreq)) :
      mcFile = findpvFile( linefreq[n], searchType=".mc" ) 
      if mcFile :
        Sdv,nch,vel,TB = getW( mcFile, [-35.,45.] )   # get integrated emission; units must be K-km/sec
      if nn > len(npList) :
        pyplot.savefig( pp, format='pdf' )
        pyplot.show()
        nn = 0
      p = pyplot.subplot(nrowpage,ncolpage,npList[nn])
      ymin = -15.
      ymax = 40.
      plotParams["ymin"] = ymin
      plotParams["ymax"] = ymax
      p.axis( [vmin, vmax, ymin, ymax], fontsize=6 )
      p.grid( True, linewidth=.01, color="0.8" )   # bigger color numbers give lighter grays!
      p.plot( vel, TB, color="r", linewidth=1. )
      p.axhline( y=0, linestyle='--', color='blue')
      p.tick_params( axis='both', which='major', labelsize=6, length=1 )
      if npList[nn]%ncols == 1 : 
        p.set_ylabel("TB (K)", fontsize=6)
      if npList[nn]-(nrows-1)*ncols > 0  : 
        p.set_xlabel("V$_{LSR}$ (km/sec)", fontsize=6)
  
    # annotate plot, shade the line
      p.text(.04, .9,  "%s" % name[n], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
      p.text(.97, .9,  "%.4f GHz" % linefreq[n], horizontalalignment='right', transform=p.transAxes, fontsize=6 )
      p.text(.04, .8,  "T$_U$ = %.0f K" % TU[n], horizontalalignment='left', transform=p.transAxes, fontsize=6 )
      p.fill_between( vel, 0., TB, color='orange', alpha=0.5, linewidth=0, \
              where=numpy.abs(vel-12.5) < 7.5 )

    # advance panel number at the very end
      nn = nn + 1
    
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

# plot centroids for several lines on one plot
# read in centroids from file Line.offsetTable
# use different symbol for each line
# color scale ranges from min to max velocity
# size ranges from min to max for each line
# when reading in files, drop scattered points immediately
# 25apr2016 - also save as RC ("Rotation Curve") object in pickle file "RCList" for
#    future rotation curve plot
# example:
# ori3.centroidPlot(["SiO_342.50","SiO_340.61","U_650.16","CO_342.65","SO_343.83","SiS_667.25"])

# I have hacked this to produce centroid figure for OriALMA;
# it is now specialized to use exactly 4 panels
# I explicitly give the positions of each subplot because I couldn't figure out
#   how to add the colorbar otherwise
# NOTE: this routine also creates pickled file RClist, which is used by ori3.figRotcurvefit()

#def centroidPlot( NameList=["SiO_342.50","SiO_340.61","CO_342.65","SO_343.83"]) :
#def centroidPlot( NameList=["SO_343.83"]) :
#def centroidPlot( NameList=["U_650.16"]) :
# def centroidPlot( NameList ) :

def centroidPlot( NameList=["SiO_342.50","SiO_340.61","CO_342.65","SO_343.83"]) :
    PA = math.pi * 142./180.  # PA of disk is 140 degrees
    cutoff = 0.04 
    symbol = [ "o", "v", "^", "s", "D", "h" ]

  # only way I could find to put colorbar on figure was to explicitly create axes
  # note that add_axes uses [left,bottom,width,height]
    figpos = [ [.1,.58,.32,.32],[.5,.58,.32,.32],[.1,.2,.32,.32],[.5,.2,.32,.32] ]
    nsym = -1
    nplot = 0

    pyplot.ioff()
    pp = PdfPages("centroids.pdf")
    fig = pyplot.figure( figsize=(8,8) )

  # diskx = 0.115 * cos(52 deg); disky = .115 * sin(52 deg)
    diskx = [ .0708, -.0708 ]
    disky = [ -.0906, .0906 ]

    for Name in NameList:
      i = Name.find('_')
      freq = float( Name[i+1:] )     
      if freq < 500. :
        fname = "/big_scr6/plambeck/350GHz/miriad/" + Name + ".offsetTable"
      else :
        fname = "/big_scr6/plambeck/650GHz/miriad/" + Name + ".offsetTable"
      try :
        fin = open( fname, "r" )
      except :
        print "couldn't open file ", fname
      print "\n\n", fname
      nsym = nsym + 1
      nplot = nplot + 1
      if nsym >= len(symbol) : nsym = 0
      v = []
      x = []
      y = []
      amp = []
      xerr = []
      yerr = []
      for line in fin :
        a = line.split()
        ampl = float( a[2] )
        xerror = float( a[4] )
        yerror = float( a[6] )
        if (ampl > 0.01) and (xerror < cutoff) and (yerror < cutoff) :
          #print line
          v.append( float(a[0]) )
          amp.append( ampl )
          x.append( float(a[3]) )
          y.append( float(a[5]) ) 
          xerr.append( xerror )
          yerr.append( yerror )
      fin.close()

    # rotate to new frame, save as RC ("Rotation Curve") object for future rotation curve plot 
      RC = {}
      RC["Name"] = Name
      RC["v"] = v
      RC["amp"] = amp

    # these are the lines I had originally; they are incorrect
      #RC["p"] = math.cos(PA)*numpy.array(x) + math.sin(PA)*numpy.array(y)
      #RC["dp"] = math.cos(PA)*numpy.array(xerr) + math.sin(PA)*numpy.array(yerr)

      RC["p"] = math.sin(PA)*numpy.array(x) + math.cos(PA)*numpy.array(y)
      RC["dp"] = math.sin(PA)*numpy.array(xerr) + math.cos(PA)*numpy.array(yerr)
      for nn in range(0,len(x)) :
        print "%8.2f  %6.3f  %6.3f  %6.3f"  % (v[nn],x[nn],y[nn],RC["p"][nn])
      fout = open( "RClist", "ab" )     # binary because I'll be pickling to it; append because I'll append
      pickle.dump( RC, fout )
      fout.close() 
      
    # size array is normalized so largest point is always the same color for each symbol
      s = numpy.array( amp )
      s = s * 200./s.max()
      p = fig.add_axes( figpos[nplot-1] )      # main plot
      #p = pyplot.subplot(3,2,nplot, aspect=1)
      p.plot( diskx, disky, linestyle="--", color='black', linewidth=4 ) 
      p.scatter( [0.], [0.], marker='*', color='black', edgecolors='black',s=300, linewidths=0.3 )
      p.axis( [.12, -.12, -.12, .12], fontsize=6 )
      p.tick_params( axis='both', which='major', labelsize=10 )
      p.errorbar( x, y, xerr=xerr, yerr=yerr, fmt='none', ecolor='black', elinewidth=.2, capsize=0.  )
      p.scatter( x, y, cmap='rainbow', c=v, marker='o', s=s, vmin=-13., vmax=25., linewidths=0.1 )
          # use edgecolor='none' to get rid of black edges
          # use s=s to make symbol areas pmoportional to flux
          # use marker=symbol[nsym] for different symbols
      p.text(.05, .91,  "%s" % Name, horizontalalignment='left', transform=p.transAxes, fontsize=12 )
      if nplot%2 == 1 :
        p.set_ylabel("$\Delta$dec (arcsec)", fontsize=10 )
      if nplot > 2 :
        p.set_xlabel("$\Delta$RA (arcsec)", fontsize=10 )

    ax2 = fig.add_axes([.3,.1,.32,.02])    # color wedge
    norm = matplotlib.colors.Normalize(vmin=-13.,vmax=25.)
    cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap='rainbow', norm=norm, orientation='horizontal' )
    cb1.set_label("LSR velocity (km/sec)", fontsize=10)
    cb1.ax.tick_params( labelsize=10 )
    pyplot.savefig( pp, format='pdf',size=10 )
    pp.close()
    pyplot.show()

# plot rotation curve from RClist data
# note: be careful! in my haste I am just putting the labels into the code!!

def RCplot( rotcurveList=["rotcurve_5Mo.dat","rotcurve_10Mo.dat"], labelList = ["5 Mo model","10 Mo model"] ) :
    color = ["red","blue","green","orange","purple","cyan","yellow"]
    RClist = []   # list of RC dictionaries
    fin = open( "RClist", "rb" )
    while (True) :
      try :
        RC = pickle.load( fin )
        RClist.append( RC )
        print RC["Name"]
      except :
        break 
    print "read in %d RC objects" % ( len(RClist) )
    fin.close()

    pyplot.ioff()
    pp = PdfPages("snippets.pdf")
    fig = pyplot.figure( figsize=(8,11) )
    p = pyplot.subplot(1,1,1)
    p.axis( [.12, -.095, -22.,32.], fontsize=12 )
    ncol = -1
    p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
    for RC in RClist :
      ncol = ncol + 1
      First = True
      if ncol > len(color) : ncol = 0
      for v,x,dx in zip( RC["v"],RC["p"],RC["dp"] ) :
        if First :
          p.plot( [x-dx,x+dx], [v,v], linewidth=6, alpha=0.5, color=color[ncol],  label=RC["Name"] )
          First = False
        else :
          p.plot( [x-dx,x+dx], [v,v], linewidth=6, alpha=0.5, color=color[ncol] )
        #p.scatter(x,v ) 

  # add theoretical rotation curves, if "rotcurve" file exists

    linestyleTable=["-","--"]
    for rotcurveFile,linestyle,label in zip(rotcurveList,linestyleTable,labelList) :    
      fin2 = open(rotcurveFile,"r")
      vr = []
      xr = []
      for line in fin2 :
        if not line.startswith("#") :
          a = line.split()
          vr.append( float(a[0]) )
          xr.append( float(a[2]) )
      fin2.close()
      p.plot( xr, vr, linewidth=1.8, color="black", label=label, linestyle=linestyle )

    p.set_xlabel("offset (arcsec)", fontsize=12)
    p.set_ylabel("LSR velocity (km/sec)", fontsize=12 )
    p.legend( loc=2, prop={'size':12}, fancybox=True, handlelength=3.2, borderaxespad=4 )
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

# compute velocity vs offset using Hirota et al. 2014 model 1
# resolve ring into 10 ringlets; compute max velocity of each ringlet;
#   for any given velocity, compute contribution from each ring,
#   then form weighted avg to get offset

# def vmodel( rin=20, rout=45, xoff=0., voff=8., vtherm=4, Mstar=5, Mdisk=0.01 ) :
def vmodel( rin=20, rout=50, xoff=0., voff=6., vtherm=4, Mstar=5, Mdisk=0.0001, rotcurveFile="rotcurve_5Mo.dat" ) :
    G = 6.672e-8    # cm3 g-1 sec-1

  # create velocity array spanning -30 to +30 km/sec, and matching amp and xmom arrays
    dv = 0.2
    vobs = numpy.arange( -30., 30.+dv, dv)
    amp = numpy.zeros( len(vobs), dtype=float  )
    xmom = numpy.zeros( len(vobs), dtype=float )

  # thermal broadening coefficients cover velocity range 0 to 2*vtherm in steps of dv
    deltav = numpy.arange( 0., 2.*vtherm+dv, dv )
    ath = numpy.exp( -2.773 * pow(deltav/vtherm, 2.))    # exp(-4 ln2 (deltaV/vtherm)^2)
    print deltav
    print ath

  # create 200-element theta array that encompasses 2pi
    dtheta = 0.005 * 2 * math.pi
    theta = numpy.arange(0.00, 1.005, 0.005) * 2.*math.pi
       # this is angle of central ray through each area element
    sintheta = numpy.sin(theta)
       # this is sin(theta) for each array element

  # step through series of rings; compute amp and v for each azimuthal element
    dr = (rout-rin)/20.
    for r in numpy.arange( rin+dr/2., rout, dr ) :
      M = 1.99e33 * (Mstar + Mdisk * (pow(r,2) - pow(rin,2))/(pow(rout,2) - pow(rin,2)))
	     # enclosed mass = central star + fraction of disk mass enclosed
      vrot = -1.e-5 * numpy.sqrt( G * M / (r*1.5e13) )   
         # Keplerian velocity in km/sec of gas at radius r
      xarr = r * sintheta        # array of x offsets
      varr = vrot * sintheta     # array of velocities
      area = r * dtheta * dr  # area per pixel for this radius

    # add contribution from each elemental vector into amp and xmom arrays
    # note: could crash if thermal broadening takes us out of vobs range
      for x1,v1 in zip(xarr,varr) :
        n = round( (v1 - vobs[0])/dv )   # find nearest vobs channel
        if n < len(vobs) :
          amp[n] = amp[n] + ath[0] * area
          xmom[n] = xmom[n] + ath[0] * area * x1
        for i in range( 1, len(ath) ) :
          if n+i < len(vobs) :
            amp[n+i] = amp[n+i] + ath[i] * area
            xmom[n+i] = xmom[n+i] + ath[i] * area * x1
          if n-i < len(vobs) :
            amp[n-i] = amp[n-i] + ath[i] * area
            xmom[n-i] = xmom[n-i] + ath[i] * area * x1
          
  # renormalize offset array
    fout = open(rotcurveFile,"w")
    fout.write("# output of ori3.vmodel\n")
    fout.write("# rin = %.2f AU\n" % rin)
    fout.write("# rout = %.2f AU\n" % rout)
    fout.write("# Mstar = %.2f Mo\n" % Mstar)
    fout.write("# Mdisk = %.2f Mo\n" % Mdisk)
    fout.write("# voff = %.3f km/sec\n" % voff)
    fout.write("# vtherm = %.2f km/sec\n" % vtherm)
    for n in range(0, len(vobs)) :
      if amp[n] > 0. :
        xmom[n] = xmom[n]/amp[n]
        fout.write( "%8.3f %8.4f %8.4f\n" % (voff+vobs[n], amp[n], xoff+(xmom[n]/415.)) )    # convert to arcsec
    fout.close()

# eqn 1,2,3 in hirota et al 2014
def hirotafit() :
  for logf in numpy.arange(1.,4.05,.05):
    fHz = pow(10.,logf) * 1.e9
    logS1 = -1.51 + 1.60*logf
    S1 = pow(10.,logS1)
    S2 = ((.082 * 2.e33)/pow( 415*3.086e18,2) * 0.1*pow(fHz/1.2e12,1.35) * \
      2.*hplanck*pow(fHz,3.)/pow(ccgs,2.) * 1./(math.exp(hplanck*fHz/(kboltzmann*100.)) -1.))/1.e-26
          # dividing by 1.e-26 to convert from cgs to milliJy
    logS2 = math.log10(S2)
    print logf, logS1, logS2, math.log10(S1+S2)

# this is what I used to make figure "rotcurvefit.pdf"
def figRotcurvefit() :
  vmodel( rin=20., rout=50, xoff=0, voff=8., vtherm=4, Mstar=5, Mdisk=0.0001, rotcurveFile="rotcurve_5Mo.dat" ) 
  vmodel( rin=20., rout=50, xoff=0., voff=8., vtherm=4, Mstar=10, Mdisk=0.0001, rotcurveFile="rotcurve_10Mo.dat" ) 
  #vmodel( rin=20., rout=50, xoff=0.01, voff=6.8, vtherm=4, Mstar=5, Mdisk=0.0001, rotcurveFile="rotcurve_5Mo.dat" ) 
  #vmodel( rin=20., rout=50, xoff=0.01, voff=6.8, vtherm=4, Mstar=10, Mdisk=0.0001, rotcurveFile="rotcurve_10Mo.dat" ) 
  RCplot( rotcurveList=["rotcurve_5Mo.dat","rotcurve_10Mo.dat"] )

# plotFlux is custom-designed to plot SrcI spectral energy distribution
# SrcIFlux.dat is nominally a csv file downloaded from Google Drive; first few columns are
#   freq,flux(mJy),fluxUncertainty,symbolsize,symbolcolor,symboltype,...
# for some maddening reason, this program seems to skip over one of the columns unless I
#   cut and paste this csv file to another file; that's what I have been doing
# 
def plotFlux( FluxFile="SrcIFlux.dat") :

  # begin by reading measured fluxes from FluxFile; expect csv format
    freq = []
    flux = []
    unc = []
    size = []
    color = []
    mtype = []
    fin = open( FluxFile, "r" )
    for line in fin :
      if not line.startswith("#") :
        a = line.split(',')
        if (len(a) > 2) :
          print a
          freq.append( float(a[0]) )
          flux.append( float(a[1]) )
          unc.append( float(a[2]) )
          size.append( 10.*float(a[3]) )
          color.append( a[4] )
          if a[5] == 't' : a[5] = '^'
          mtype.append( a[5] )
    fin.close()
    print color
      # this is here because sometimes a column is skipped?!
    print mtype

  # begin log-log plot
    pyplot.ioff()
    pp = PdfPages("Flx.pdf")
    pyplot.figure( figsize=(8,8) )
    ax = pyplot.subplot(1,1,1)
    # ax.axis( [4, 1000, .1, 20000.], fontsize=12, linewidth=1.1 )
    ax.axis( [4, 1000, .1, 20000.], fontsize=20, linewidth=1.1 )
    ax.set_xscale( "log" )
    ax.set_yscale( "log" )
    # ax.set_xlabel("freq (GHz)", fontsize=12)
    ax.set_xlabel("freq (GHz)", fontsize=16)
    # ax.set_ylabel("flux density (mJy)", fontsize=12)
    ax.set_ylabel("flux density (mJy)", fontsize=16)
    ax.tick_params(length=8, which='major', labelsize=16)
    ax.tick_params(length=6, which='minor', labelsize=16)
    pyplot.grid(True)

  # plot nu^2 curve; show as dashed line below 30 GHz
  # nu^2 curve passes through 340 mJy at 229 GHz
    ax.plot([1,30.],[.00629,5.6635],linestyle='--',color='black',linewidth=0.8)
    ax.plot([30.,1000.],[5.6635,6293.],linestyle='-',color='black',linewidth=1.2)

  # plot Hminus curve from file Hminus.dat, if it can be found
    hfreq = []
    hflux = []
    try :
      fin = open("Hminus.dat")
      for line in fin :
        if not line.startswith('#') :
          a = line.split()
          hfreq.append( float(a[0]) )
          hflux.append( float(a[3]) )
      fin.close()
      ax.plot( hfreq, hflux, linestyle='--', color='red', linewidth=1.2)
    except:
      pass

  # now plot the data points one by one with errorbar, which does not seem to accept lists of styles, colors, sizes
    for frq,flx,un,siz,col,mtyp in zip(freq,flux,unc,size,color,mtype) :
      pyplot.errorbar( frq, flx, yerr=un, fmt=mtyp, capsize=0, markersize=0.8*siz, color=col, elinewidth=2, ecolor='black' )
    # ax.plot([4.74,7.34],[.3904,.88],linestyle=':',color='green')
      # this shows spectral index derived by Forbrich et al

  # label the dust and free-free curves
    #ax.text( 0.07, .9, "Orion SrcI", transform=ax.transAxes, \
    #    horizontalalignment='left', fontsize=18, rotation='horizontal' )
    ax.text( 0.66, .31, "free-free", transform=ax.transAxes, \
        horizontalalignment='center', fontsize=16, color="red", rotation='horizontal' )
    ax.text( 0.66, .56, "dust", transform=ax.transAxes, \
        horizontalalignment='center', fontsize=16, color="black", rotation=41 )

  # brute force way to add legends; symbol size scheme is different, but this looks OK
  # another problem is that a line is drawn through symbols on the legend, but it is so short that it mostly doesn't show
    ax.plot( .01, .01, 'bs', ms=8, label="Plambeck 2013")
    ax.plot( .01, .01, marker='d', ms=8, color='cyan', label="Zapata 2004")
    ax.plot( .01, .01, marker='o', color='yellow', ms=9, label="Beuther 2004,2006")
    ax.plot( .01, .01, marker='d', color='orange', ms=10, label="Hirota 2014,2015")
    ax.plot( .01, .01, marker='^', color='magenta', ms=8, label="Rivilla 2015")
    ax.plot( .01, .01, 'gD', ms=8, label="Forbrich 2016")
    ax.plot( .01, .01, 'rs', ms=10, label="this paper")
    ax.legend( numpoints=1, loc=2, prop={'size':10}, markerscale=1, fancybox=True, handlelength=1, borderaxespad=3 )
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

# deconvolve source size from measured size and beamwidth
#def decon( beamMaj, beamMin, beamAngle, srcMaj, srcMin, srcAngle, pa=140) :
  # first find beam diameter at position angle pa
#    thetaBm = o

# find diameter of ellipse at angle pa; use eqn 10.23 in math handbook
def dEllipse( major, minor, majDeg, paDeg ) : 
    theta = math.radians( majDeg - paDeg )
    rsq = pow(major*minor,2)/( pow(major*math.sin(theta), 2) + pow(minor*math.cos(theta), 2) )
    print (majDeg - paDeg), math.sqrt(rsq)

# generate olay file from ForbrichCatalog.txt
def Forbrich() :
    fin = open( "/o/plambeck/OriALMA/Refs/ForbrichCatalog.txt", "r" )
    fout = open( "forbrich.olay", "w" )
    for line in fin :
      if not line.startswith("#") :
        a = line.split()
        fout.write("sym hms dms %s no %s %s %s %s %s %s 2 1\n" % (a[0], a[2], a[3], a[4], a[6], a[7], a[8]))
    fin.close()
    fout.close()
      
# check our solution to eqn 3 for the referee
def check3( ) :
    for tauff in numpy.arange(1.,5.,.5) :
      tauL = 5.7*tauff
      print "%8.1f  %8.4f" % ( tauff, (1.-math.exp(-1.*(tauff+tauL)))/(1.-math.exp(-1.*tauff)))

# make list of cuts through SiS map - for velplot and for cgdisp olay
# note: x coordinate into velplot has reversed sign (negative is to the left)
#   relative to convention in cgdisp, for arcsec offsets
def cuts( olayFile="cutsolay", velplotfile="pvcmds" ) :
  pacut = 140
  pastep = pacut - 90.
  fout = open( olayFile, "w" )
  n = 0
  #for dist in numpy.arange(-2.8,3.0,.4) :
  for dist in numpy.arange(-1.4,1.5,.2) :
    n = n+1
    dx = dist * math.sin( math.radians(pastep) )
    dy = dist * math.cos( math.radians(pastep) )
    fout.write("vector arcsec arcsec %d yes %.3f %.3f 100 %d\n" % (n,dx,dy,pacut))
    fout.write("vector arcsec arcsec %d no %.3f %.3f 100 %d 0\n" % (n,dx,dy,pacut-180.))
    print "%d, %.3f, %.3f, %d" % (n,-dx,dy,pacut)
  fout.close()

def plotcuts( smin=-.02, smax=.08, nx=2, ny=2, nslice1=1, nslice2=15 ) :
    pmin = -2.
    pmax = 2.
    vmin = -20
    vmax = 30 
    npanel = 0 
    pyplot.ioff()
    fig = pyplot.figure( figsize=(11,8) )
    pp = PdfPages("pvcuts.pdf")
    for nslice in range(nslice1,nslice2+2) :
      npanel = npanel + 1
      if npanel > nx*ny :
        pyplot.savefig( pp, format='pdf' )
        pyplot.show()
        pyplot.figure( figsize=(11,8) )
        npanel = 1
      print "npanel = %d" % npanel 
      p = pyplot.subplot(nx,ny,npanel)
      if nslice <= nslice2 :
        [vel, pos, flx] = readpv( "pv%d.mp" % nslice, pmin=pmin, pmax=pmax, vmin=vmin, vmax=vmax )
        nv = len( numpy.unique( vel ) )
        np = len( numpy.unique( pos ) )
        p.axis( [vel[0], vel[-1], pos[0], pos[-1]] )
        p.tick_params( which='major', labelsize=6)
        p.tick_params( which='minor', labelsize=6)
        imgplot = p.imshow(numpy.reshape(flx, (np,nv) ), origin='lower', aspect='auto', \
            extent=[vel[0],vel[-1],pos[0],pos[-1]], vmin=smin, vmax=smax) 
        p.grid( True, linewidth=0.1, color='0.1' )   # color=0.1 is a light gray
        p.text(.04, .88,  "slice %d" % nslice, horizontalalignment='left', transform=p.transAxes, fontsize=8, color='white' )
      else :
        pyplot.axis('off')
        cbar = pyplot.colorbar( imgplot )
        cbar.ax.tick_params(labelsize=6)

  # extract original file name from velplot history item
    label = "no label"
    pr = subprocess.Popen( ( shlex.split("head pv%d.mp/history" % nslice2 ) ), \
      stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
    result = pr.communicate()[0]
    n1 = str.find( result, "in=" )
    n2 = str.find( result[n1:], "VELPLOT" )
    if n1 >= 0 :
      label = result[n1+3:n1+n2-1]
    pyplot.suptitle( label )
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

def archiveSrch( infile ) :
    fin = open(infile,"r")
    projList = []
    RACol = 0
    DecCol = 0
    ResolutionCol = 0
    ArrayCol = 0
    MosaicCol = 0
    FreqSupportCol = 0
    for line in fin:
      if not line.startswith("#") :
        a = line.split("\t")

      # figure out which column is which based on 1st line
        if "Project code" in a[0] :
          CodeCol = 0
          for n in range(1,len(a)) :
            if "RA" in a[n] :
              RACol = n
            if "Dec" in a[n] :
              DecCol = n
            if "Release date" in a[n] :
              ReleaseCol = n
            if "Spatial resolution" in a[n] :
              ResolutionCol = n
            if "Array" in a[n] :
              ArrayCol = n
            if "Mosaic" in a[n] :
              MosaicCol = n
            if "Largest angular scale" in a[n] :
              largestAngScaleCol = n
            if "Frequency support" in a[n] :
              FreqSupportCol = n
            if "PI name" in a[n] :
              PIcol = n
            if "Field of view" in a[n] :
              FovCol = n
      # create dictionary holding values; append to list of dictionaries
        else :
          if len( a[ReleaseCol] ) < 2 :     # in a few cases, field is blank
            a[ReleaseCol] = "3000-01-01"
          proj = { 'code' : a[CodeCol], 
                   'spw'  : spw( a[FreqSupportCol] ),
                   'res'  : float( a[ResolutionCol] ),
                   'LAS'  : float( a[largestAngScaleCol] ),
                   'release' : a[ReleaseCol],
                   'PI'   : a[PIcol],
                   'FOV'  : a[FovCol],
                   'RA'   : a[RACol],
                   'DEC'  : a[DecCol]
                 }
          projList.append( proj )
            # proj('spw') = spw
            # spw[y1][f1,f2,dv] = fstartGHz,fstopGHz,dvkms for spectral window y1
    fin.close()
    return projList
          
# spw decodes the Frequency support column
# returns list of spectral windows; each one [f1GHz,f2GHz,dvkms]
def spw( s ) :
    a = s.split("[")
    spw = []
    nn = 0
    for n in range(0,len(a)) :
      if len(a[n]) > 0 :
        b = a[n].split(",")
        n1 = string.find( b[0], ".." )
        n2 = string.find( b[0], "GHz" )
        spwstart = float( b[0][0:n1] )
        spwstop = float( b[0][n1+2:n2] )
        n3 = string.find( b[1], "kHz" )
        dfGHz = float( b[1][0:n3] )/1.e6   # convert channel width from kHz to GHz
        dv = 2.998e5 * dfGHz/((spwstart+spwstop)/2.)
        spw.append( [spwstart,spwstop,dv] )
    return spw

# convenience function to see if spectral window spw lies within frqBand
def inrange( spw, frqBand ) :
    if (spw[0] < frqBand[0]) and (spw[1] < frqBand[0]) :    # entire spw is below freq range
      return False
    elif (spw[0] > frqBand[1]) and (spw[1] > frqBand[1]) :    # entire spw is above freq range
      return False
    else :
      return True	  # some piece of spw must lie inside frqBand

# convenience function to turn RA or DEC ascii triplet into a decimal
def decimal( triplet ) :
   a = triplet.split(":")
   if float(a[0]) > 0. :
     return float(a[0]) + float(a[1])/60. + float(a[2])/3600.
   else :
     return float(a[0]) - float(a[1])/60. - float(a[2])/3600.
 
# convenience function to see if radec of target lies within FOV of observation
def inbeam( radec, proj ) :
    RAtarg = 15.*decimal( radec[0] )
    DECtarg = decimal( radec[1] )
    RAobs = float( proj['RA'] )     # it appears that archive gives decimal degrees for RA
    DECobs = float( proj['DEC'] )
    deltaRA = (RAtarg - RAobs)/math.cos(math.radians(DECtarg))
    deltaDEC = DECtarg - DECobs
    offset = 3600*math.sqrt( pow(deltaRA,2) + pow(deltaDEC,2) ) #/ float(proj['FOV'] )
       # offset, expressed as fraction of field of view; e.g. 0.5 means at half pwr point
       # mosaics are a problem; not sure what to do with these 
    print "offset = %.2f, FOV = %.2f" % (offset, float(proj['FOV']))
    if offset < 0.6*float(proj['FOV']) :
       return True
    else :
       return False
        
# generate plot of frequency coverage for data in ALMA archive
# infile = tab separated output from search query
# frqBand = [f1GHz,f2GHz] = freq range of interest (plot only datasets within this range)
#
recombList = [85.688,92.034,99.023,106.737,124.747,135.286,147.047,160.211,174.995,191.657,210.502, \
  231,901,256.302,284.250,316.415,353.623,396.901,447.540,507.175,577.896,662.404]

def archivePlot( infile, frqBand, radec=["05:35:14.109","-05:22:22.73"], frqMarkerList=recombList ) :
    projList = archiveSrch( infile ) 
    today = datetime.datetime.today()
    deltavert=.03
    rheight = .02
    pyplot.ioff()
    pp = PdfPages("projects.pdf")
    fig = pyplot.figure( figsize=(11,8) )
    p = fig.add_axes( [.05,.08,.6,.9] ) 
    pyplot.axis( [frqBand[0],frqBand[1],0.,1.] )
    vert = 0.                # allowed range is 0-110
    for proj in projList :
      newproj = True
      for spw in proj['spw'] :
        if inrange( spw,frqBand ) and inbeam( radec, proj ) :
          if newproj:         # allow extra 1-unit gap between projects
            vert = vert+deltavert
            if vert > (1.-deltavert/2.) :
	          print "WARNING: plot is full"
            newproj = False 
            codeColor = "black"
            print  proj['release'] 
            if dateutil.parser.parse( proj['release'] ) > today :
              codeColor = "red"    # data not yet available
            pyplot.text(1.01, vert+rheight/2., proj['code'], horizontalalignment='left', \
              verticalalignment='center',transform=p.transAxes, color=codeColor, fontsize=10 )
            pyplot.text(1.26, vert+rheight/2., "%5.2f" % proj['res'], horizontalalignment='right', \
              verticalalignment='center',transform=p.transAxes, fontsize=10)
            pyplot.text(1.34, vert+rheight/2., "%5.1f" % proj['LAS'], horizontalalignment='right', \
              verticalalignment='center',transform=p.transAxes, fontsize=10)
            ncomma = string.find( proj['PI'], "," )
            pyplot.text(1.37, vert+rheight/2., proj['PI'][0:ncomma], horizontalalignment='left', \
              verticalalignment='center',transform=p.transAxes, fontsize=10)
          vcolor = "gray"
          if spw[2] < 3. :
            vcolor = "blue"
          if spw[2] < 1. :
            vcolor = "red"
          rect = Rectangle( (spw[0],vert), (spw[1]-spw[0]), rheight, facecolor=vcolor, alpha=0.4,\
            edgecolor=vcolor, linewidth=1 ) 
          pyplot.axhline( y=vert+rheight/2., color='gray', linestyle='dotted')
          p.add_patch( rect ) 
    pyplot.yticks([])
    if len(frqMarkerList) > 0 :
      for frq in frqMarkerList :
        pyplot.axvline( x=frq, color="red", linestyle="dashed")
    p.grid(True)
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

#clumpList = [ { "name": "SrcI",    "region":"arcsec,box(0.50,-0.50,-0.50,0.50)", "vlsr" : 5. }, \
#              { "name": "HC1",     "region":"arcsec,box(1.49,-0.73,0.49,0.27)",  "vlsr" : 5. }, \
#              { "name": "HC2",     "region":"arcsec,box(0.59,-2.22,-0.41,-1.22)", "vlsr" : 5. }, \
#              { "name": "HC3",     "region":"arcsec,box(-0.76,-4.93,-1.76,-3.93)", "vlsr" : 5. }, \
#              { "name": "BN",      "region":"arcsec,box(-5.57,7.35,-6.57,8.35)", "vlsr" : 21. }, \
#              { "name": "IRc7",    "region":"arcsec,box(-3.34,-1.34,-4.34,-0.34)", "vlsr" : 10. }, \
#              { "name": "W1",      "region":"arcsec,box(-4.96,2.07,-5.96,3.07)", "vlsr" : 10. }, \
#              { "name": "SW",      "region":"arcsec,box(-5.71,-6.93,-6.71,-5.93)", "vlsr" : 10. }, \
#              { "name": "N",       "region":"arcsec,box(6.29,12.57,5.29,13.57)", "vlsr" : 10. }, \ 
#              { "name": "N2",      "region":"arcsec,box(1.7,.7,.7,1.7)", "vlsr" : 10. } ]
#  
#              { "name":"SrcIdisk", "region":"arcsec,polygon(.14,-.2,-.03,.02,.01,.07,.22,-.14)",  "vlsr":5.}   use for highres

# added 9 nov 2017; trying to find region with simpler spectrum for Suchi to fit
clumpList = [ { "name": "SW",   "region" : "arcsec,box(-5.71,-6.93,-6.71,-5.93)", "vlsr" : 10. } ] 

# use imlist options=stat to get statistics of region; return [plane, max, min, rms] arrays
def dumpstats( infile, region, tmpfile="junk" ) :
    plane = []
    smax = []
    smin = []
    srms = [] 
    p= subprocess.Popen( ( shlex.split("imlist in=%s region=%s options=stat log=%s" % \
      (infile,region,tmpfile) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    time.sleep(1)
    result = p.communicate()[0]
    print result
    if "Fatal Error" in result :
      print " --- fatal --- "
      return
    n = 0
    fin = open( tmpfile, "r" )
    for line in fin :
      n = n+1
      if n>12 :
        a = line.split()
        plane.append( int(a[0]) )
        smax.append( float(a[3]) )
        smin.append( float(a[4]) )
        srms.append( float(a[6]) )
    return plane, smax, smin, srms
    fin.close()

# mosspectra is designed to write out spectra for Suchi's project
def mosspectra( ) :
    noisebox = 'arcsec,box(10)'
    for clump in clumpList:
      for spw in ["spw0", "spw1", "spw2", "spw3"] :
        infile = "/alma_scr/plambeck/220GHz_mosaic/miriad/mos.%s.cm" % spw
        hann = 1
        # outfile = "/o/plambeck/OriALMA/220Spectra/%s.%s.txt2" % (clump["name"],spw)
        outfile = "/alma_scr/plambeck/XCLASS/%s.%s.dat" % (clump["name"],spw)
        print infile
        [chan, freq, flux ] = getspec( infile, region=clump["region"], \
           vsource=clump["vlsr"], hann=hann )
        [plane, smax, smin, srms ] = dumpstats( infile, region=noisebox )
        print "spectra contains %d chans; stats contains %d planes" % (len(chan),len(plane))
        if len(plane) != len(chan) :
          print "FATAL ERROR"
          return
        fout = open( outfile, "w" )
        fout.write("# spectrum created with ori3.mosspectra\n")
        fout.write("# infile: %s\n" % infile)
        fout.write("# region: %s\n" % clump["region"] )
        fout.write("# vsource: %.2f\n" % clump["vlsr"] )
        fout.write("# hann: %d\n" % hann)
        fout.write("# noisebox: %s\n" % noisebox )
        fout.write("#  freq_MHz       flux     map_max     map_min    map_rms   plane\n")
        for [chn,frq,flx] in zip( chan, freq, flux ) :
          n = int(chn)-1 
          fout.write("%13.5f  %9.7f  %9.7f  %9.7f  %9.7f  %5d\n" % (1000.*frq,flx,smax[n],smin[n],srms[n],chn))
        fout.close() 

# short routine to read absolute positions from olay file, generate regions
def clumpBoxes( infile="xx.olay", box=1. ) :
    fin = open(infile, "r" )
    for line in fin :
      a = line.split()
      raSec = float(a[7])
      decSec = float(a[10])
      draArcsec = 15.*(raSec - 14.514)          # RA offset from SrcI
      ddecArcsec = 30.575 - decSec     # negative ddec for larger DEC because decs are negative
      print "(%.2f,%.2f,%.2f,%.2f)" % (draArcsec+box/2.,ddecArcsec-box/2.,draArcsec-box/2.,ddecArcsec+box/2.)

# generate vector field as overlay file (to compare SiO and 29SiO pol directions on one plot, 5/10/18)
def vecs( poli, pa, region, olay, delta=.03, scale=5., color=0, writemode="w" ) :
  
    p = subprocess.Popen( ( shlex.split("imtab in=%s log=imtablog_poli format=(3F12.5) region=%s" % (poli,region) ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
    result1 = p.communicate()[0]
    a1 = result1.split()
    p = subprocess.Popen( ( shlex.split("imtab in=%s log=imtablog_pa format=(3F12.5) region=%s" % (pa,region) ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
    result2 = p.communicate()[0]
    a2 = result2.split()
    print a2[6],a2[7],a2[8]
    if a2[7] == "Fatal" :
      return
  
  # open the logs, fill out lists, then convert to array
    fin_poli = open("imtablog_poli","r")
    fin_pa = open("imtablog_pa","r")
    fout = open(olay,writemode)
    fout.write("# overlay file created by ori3.vecs\n")
    fout.write("# poli file: %s\n" % poli)
    fout.write("# pa file: %s\n" % pa)
    fout.write("# region: %s\n" % region)
    fout.write("# delta = %.3f (gap between vectors, arcsec)\n" % delta)
    fout.write("# scale = %.3f (scale - negative indicates all one length)\n" % scale)
    fout.write("# each segment is created with 2 vectors starting at sampling point, 180 degrees apart\n")
    fout.write("color %d\n" % color)
    for line1,line2 in zip(fin_poli,fin_pa) :
      a = line1.split()
      b = line2.split()
      if (a[0] != b[0]) or (a[1] != b[1]) :
        print "error - positions don't match"
      else :
      # testing ( x % delta < 1.e-5 ) sometimes failed because it returned ~delta rather than 0
        fracRA,wholeRA = math.modf( float(a[0])/delta )
        fracDEC,wholeDEc = math.modf( float(a[1])/delta )
        if abs(fracRA) < 0.1*delta and abs(fracDEC) < 0.1*delta :
          #print a[0],a[1],a[2],b[2]
          if scale > 0.:
            length = scale*float(a[2])
          else :
            length = -1.*scale
        # cgdisp documentation is wrong; x,y is not center of vector, it is starting point;
        # therefore, write out 2 half-length vectors 180 degrees apart
          fout.write("vector arcsec arcsec vec no %s %s %8.4f %s 0\n" % \
            (a[0],a[1],length/2.,float(b[2])) )
          fout.write("vector arcsec arcsec vec no %s %s %8.4f %s 0\n" % \
            (a[0],a[1],length/2.,float(b[2])-180.) )
    fin_poli.close()
    fin_pa.close()
    fout.close()

# hann designed to read stacked spectra from Adam
# writes out files with format compatible with dumpspec (chan, freqGHz, amp are 1st 3 items)
def hann( infile, outfile, hann=15 ) :
    fGHz,amp = numpy.loadtxt( infile, unpack=True )
    h = numpy.hanning( hann+2 )	# for odd n, numpy.hanning returns a window function with 0 as first and last values
    sm = h[1:-1]/numpy.sum(h)
      # normalize the array
    amp2 = numpy.convolve( amp, sm, mode="same")
    fout = open( outfile, "w" )
    fout.write("# created with ori3.hann\n")
    fout.write("# infile: %s\n" % infile)
    fout.write("# hann: %d\n" % hann)
    fout.write("# vsource: 5.0\n")    # this is only for SrcI spectra
    fout.write("# col 1 = (fake) channel, 2 = freq, 3 = smoothed spectrum, 4 = original spectrum\n")
  # write out in format that can be used by readspec
    n = 0
    for f,a,a2 in zip(fGHz,amp,amp2) :
      n = n + 1
      fout.write("%5d  %10.6f %10.5f %10.5f\n" % (n,f,a2,a))
    fout.close()

# creates histogram of position angles, does Gaussian fit to get sigma
# designed for Davis-Chandra-Fermi method, 3/28/19
def PAhisto( imageFile, region="arcsec,box(1.2)" ) :
    p = subprocess.Popen( ( shlex.split("imtab in=%s log=imtablog region=%s format=(3F12.5)" % (imageFile,region) ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
    result = p.communicate()[0]
    x,y,pa = numpy.loadtxt("imtablog", unpack=True )
    print "found %d points" % len(pa)

  # adjust points affected by 180 degree ambiguities
    nadjust = 0
    for i in range(0,len(pa)) :
      if pa[i] < -145. :
        #print "%12.5f %12.5f %12.5f" % ( x[i],y[i],pa[i] )
        pa[i] = pa[i] + 180.
        nadjust = nadjust + 1
      if pa[i] > 145. :
        #print "%12.5f %12.5f %12.5f" % ( x[i],y[i],pa[i] )
        pa[i] = pa[i] - 180.
        nadjust = nadjust + 1
    print "corrected 180 degree shifts on %d/%d points" % (nadjust,len(pa))

  # make the plot
    pyplot.ioff()
    pp = PdfPages("histo.pdf")
    fig = pyplot.figure( figsize=(11,7) )
    n, bins, patches = pyplot.hist( pa, bins=41, range=[-20.5,20.5], alpha=0.3 )
    mu,sigma = norm.fit(pa)
    print mu,sigma
    print "sigma = %.2f degrees" % (math.sqrt(numpy.var(pa)))
    sigma = 3.
    mu = 0.
    xpdf = numpy.arange(-22.,22.01,.01)
    ypdf = 2800.*mlab.normpdf( xpdf, mu, sigma)
    pyplot.plot(xpdf,ypdf,"r-",linewidth=3)
    pyplot.tick_params(labelsize=18)
    pyplot.xlim([-22,22])
    pyplot.xlabel("$\Delta$PA (deg)", fontsize=18)
    pyplot.ylabel("N", fontsize=18)
    #pyplot.grid(True)
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

# generate data for a plot of StokesV vs poli
# use imtab to read in StokesI, poli, StokesV data on a channel-by-channel basis
# poli is masked, so will have fewer points
# find V,I pixels corresponding to each poli pixel; write [chan, velocity, x, y, I, poli, V] to outfile
# ideally this routine would check to make sure pixels, map center, etc is identical, but this version
#    does not do these checks
#
def VV( Imap, polimap, Vmap, region, chanrange ):
    
  # step 1: use imlist to retrieve velocity and freq information from the header
    p= subprocess.Popen( ( shlex.split('imlist in=%s' % Imap) ), \
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

  # now step through given channel range
    fout = open("VV","w")
    for ichan in range( chanrange[0],chanrange[1]+1) :
      vchan = v1 + (ichan-1)*dv
      p = subprocess.Popen( ( shlex.split("imtab in=%s region=%s(%d) log=imtablog format=(3F12.5)" % 
         (polimap,region,ichan) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
      result = p.communicate()[0]
      try :
        xpoli,ypoli,poli = numpy.loadtxt("imtablog", unpack=True )
        npts = len(poli)
      except :
        npts = 0
      print "%s channel %d has %d points" % (polimap,ichan,npts)

      if npts > 0 :
        p = subprocess.Popen( ( shlex.split("imtab in=%s region=%s(%d) log=imtablog format=(3F12.5)" % 
           (Imap,region,ichan) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
        result = p.communicate()[0]
        xI,yI,I = numpy.loadtxt("imtablog", unpack=True )

        p = subprocess.Popen( ( shlex.split("imtab in=%s region=%s(%d) log=imtablog format=(3F12.5)" % 
           (Vmap,region,ichan) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
        result = p.communicate()[0]
        xV,yV,V = numpy.loadtxt("imtablog", unpack=True )
 
      # the following code requires that pixel order dumped by imtab always is the same
        nI = 0
        Imax = 0.
        VImax = 0.
        VPmax = 0.
        PImax = 0.
        for n in range(0,npts) :
          while (xI[nI] != xpoli[n]) or (yI[nI] != ypoli[n]) :
            nI = nI + 1
          if (xI[nI] == xpoli[n]) and (yI[nI] == ypoli[n]) and (xV[nI] == xpoli[n]) and (yV[nI] == ypoli[n]) :
            if I[nI] > Imax :
              Imax = I[nI]
            VIratio = abs(V[nI])/I[nI]
            VPratio = abs(V[nI])/poli[n]
            PIratio = poli[n]/I[nI]
            if VIratio > VImax :
              VImax = VIratio
            if VPratio > VPmax :
              VPmax = VPratio
            if PIratio > PImax :
              PImax = PIratio 
          else :
            print "FATAL - I and V arrays don't match"
            break

      fout.write("%5d  %7.3f  %9.5f  %7.5f  %7.5f  %7.5f\n" % \
         (ichan, vchan, Imax, PImax, VImax, VPmax))
    fout.close()     

# helper function to read intensity vs velocity from imspec output
# this routine trims off leading and trailing text, then reads spectrum with numpy.loadtxt
def readImspecFile( infile ) :
    fin = open(infile, "r")
    read = False
    sp = []
    for line in fin :
      if read :
        if (len(line) < 2) :     # spectrum ends with blank line
           read = False
        else :
          sp.append(line)
      elif ("Relativist" in line ):
        read = True
    chan, vlsr, flux, dummy1, npts = numpy.loadtxt( sp, unpack=True )
    return chan,vlsr,flux
 
# reads I and V spectra produced by imspec
# computes derivative of I
# writes out new file containing [chan,vlsr,I,Ideriv,V]
def ZeemanTest( Ifile, Vfile, outfile ):
    chanI, vlsrI, fluxI = readImspecFile( Ifile ) 
    deriv = numpy.gradient( fluxI )
    chanV, vlsrV, fluxV = readImspecFile( Vfile ) 
    fout = open( outfile, "w")
    for n,v,I,Id,V in zip(chanI,vlsrI,fluxI,deriv,fluxV) :
      fout.write("%8d %8.2f %12.3e %12.3e %12.3e\n" % (n, v, I, Id, V) )
    fout.close()

def plt( x, y, z) :
    nx = len(numpy.unique(x))
    ny = len(numpy.unique(y))
    fig,ax = pyplot.subplots(nrows=2, ncols=2 )
    im = ax[0,0].imshow(numpy.reshape(z[0], (nx,ny)), origin="lower")
    ax[0,0].set_aspect('equal',adjustable='box')
    im = ax[0,1].imshow(numpy.reshape(z[1], (nx,ny)), origin="lower")
    ax[0,1].set_aspect('equal',adjustable='box' )
    im = ax[1,0].imshow(numpy.reshape(z[2], (nx,ny)), origin="lower")
    ax[1,0].set_aspect('equal',adjustable='box')
    im = ax[1,1].imshow(numpy.reshape(z[3], (nx,ny)), origin="lower")
    ax[1,1].set_aspect('equal',adjustable='box')
    pyplot.subplots_adjust( wspace=.1, hspace=.1)
    pyplot.show()
    
# compute polarization statistics (P/I, V/I, V/P) vs velocity for SiO masers
#def polStats( IQUVmapList, pcutoff, icutoff, region='arcsec,box(1.2)', chanrange=[40,100] ):
def polStats( IQUVmapList, pcutoff, icutoff, region='arcsec,box(1.2)', chanrange=[40,105] ):
    
  # use imlist to retrieve velocity information from the header
    Imap = IQUVmapList[0]
    print Imap
    p= subprocess.Popen( ( shlex.split('imlist in=%s' % Imap) ), \
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

  # write header info into each of the 6 output files
    for outfile in ["Istat.dat", "fracPstat.dat", "fracP1stat.dat", "fracVstat.dat", "fracV1stat.dat", "VoverPstat.dat"] :
      fout = open(outfile,"w")
      fout.write("# %s\n" % outfile)
      fout.write("# created by ori3.polStats\n")
      fout.write("# input files: %s\n" % IQUVmapList[0] )
      fout.write("#              %s\n" % IQUVmapList[1] )
      fout.write("#              %s\n" % IQUVmapList[2] )
      fout.write("#              %s\n" % IQUVmapList[3] )
      fout.write("# region: %s\n" % region )
      fout.write("# icutoff: %.5f\n" % icutoff )
      fout.write("# pcutoff: %.5f\n" % pcutoff )
      fout.write("#\n#  chan  vel   npixels   5_pct 25_pct 50_pct 75_pct 95_pct\n")
      fout.close()

  # read I,Q,U,V data for each channel range

    ifl = []
    fp = []
    fv = []
    vp = []

    for ichan in range( chanrange[0],chanrange[1]+1) :
      vchan = v1 + (ichan-1)*dv
      z = []
      for nmap in range( 0,4 ):
        p = subprocess.Popen( ( shlex.split("imtab in=%s region=%s(%d) log=imtablog format=(3F12.5)" % 
           (IQUVmapList[nmap],region,ichan) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
        result = p.communicate()[0]
        x,y,flx = numpy.loadtxt("imtablog", unpack=True )
        z.append(flx)

    # dump percentiles of I
      iflagged = numpy.ma.masked_where( z[0] < icutoff, z[0] )
      ifl = ifl + numpy.ma.compressed( iflagged ).tolist() 
      pdump( ichan, vchan, iflagged, "Istat.dat" )

    # compute fractional linear pol flagged where I < icutoff
      poli = numpy.sqrt(numpy.power(z[1],2) + numpy.power(z[2],2)) 
      poli = numpy.ma.masked_where( z[0] < icutoff, poli)
      polm = (poli/z[0])
      fp = fp + numpy.ma.compressed(polm).tolist() 
      pdump( ichan, vchan, polm, "fracPstat.dat" )

    # fractional linear pol flagged where I < icutoff AND P < pcutoff
      polm = numpy.ma.masked_where( poli < pcutoff, polm )
      pdump( ichan, vchan, polm, "fracP1stat.dat" )

    # fractional circular pol, flagged where I < Icutoff
      fracv = numpy.ma.masked_where( z[0] < icutoff, numpy.absolute(z[3])/z[0] )
      fv = fv + numpy.ma.compressed( fracv ).tolist()
      pdump( ichan, vchan, fracv, "fracVstat.dat" )

    # fractional circular pol, flagged where I < icutoff AND P < pcutoff
      fracv1 = numpy.ma.masked_where( poli < pcutoff, fracv )
      pdump( ichan, vchan, fracv1, "fracV1stat.dat" )

    # circular/linear pol, flagged where P < pcutoff, I < icutoff
      voverp = numpy.absolute(z[3])/poli
      voverp = numpy.ma.masked_where( z[0] < icutoff, voverp )
      vp = vp + numpy.ma.compressed( voverp ).tolist()
      voverp = numpy.ma.masked_where( poli < pcutoff, voverp ) 
      pdump( ichan, vchan, voverp, "VoverPstat.dat" )

     # plt( x, y, [z[0],polm,fracv,voverp] )
  #    print len(ifl), len(fp), len(fv), len(vp)

  # plot fracp, fracv vs I
  #  fig,ax = pyplot.subplots( nrows=2 )
  #  ax[0].scatter( ifl, fp, alpha=0.1 )
  #  pyplot.title("P/I vs I")
  #  ax[1].scatter( ifl, fv, alpha=0.1 )
  #  pyplot.title("V/I vs I")
  #  pyplot.show()
    

def pdump( ichan, vchan, arr, outfile ) :
    fout = open(outfile,"a")
    arr1 = numpy.ma.compressed( arr )      # ditch the masked values
    nl = len(arr1)
    print "%3d  %s  %d points" % (ichan, outfile, nl)
    if nl > 10 :
      a = numpy.percentile( arr1, [5,25,50,75,95] )
      fout.write("%5d  %7.2f  %5d  %7.4f %7.4f %7.4f %7.4f %7.4f\n" % \
        (ichan, vchan, nl, a[0],a[1],a[2],a[3],a[4]) )
    else :
      print "skipping percentiles for %s" % outfile
    fout.close()

# helper routine for polspec; use imspec to dump 1 spectrum from infile
def dump1( infile, region, tmpfile="junk") :
    print "reading %s" % infile
    chan = []
    vlsr = []
    y = []
    ngood = 0
    nflagged = 0
    p= subprocess.Popen( ( shlex.split("imspec in=%s region=%s options=list,eformat,noheader log=%s" % \
      (infile,region,tmpfile) )), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    time.sleep(1)
    result = p.communicate()[0]
    if "Fatal Error" in result :
      print " --- fatal --- "
      return
  # can't use numpy.loadtxt because sometimes points are masked
    fin = open(tmpfile, "r")
    for line in fin :
      if "All points masked for plane" in line :
        nflagged = nflagged + 1
        chan.append(int(line[27:])) 
        vlsr.append(-10000.)
        y.append(0.)
      else :
        ngood = ngood + 1
        a = line.split()
        chan.append(int(a[0]))
        vlsr.append(float(a[1]) )
        y.append(float(a[2]))
    fin.close()
    print "%d good values, %d flagged values" % (ngood,nflagged)
    return chan,vlsr,y

# helper routine for polspec to take out 180 degree flips in PA
# midpt is midpt of range; newpa constrained to be in the range midpt-90 to midpt+90
def paflip( pa, midpt=80 ) :
    panew = []
    for angle in pa :
      while angle > midpt+90.:
        angle = angle - 180.
      while angle < midpt-90. :
        angle = angle + 180.
      panew.append(angle)
    return panew

# dump out pickle file containing [chan, vel, I, poli, pa, paerr, V] table for a particular region
def polspec( rootfile, region, outfile, tmpfile="junk" ):
    chan,vlsr,I = dump1( rootfile + ".I.cm", region )
    dummy,dummy,poli = dump1( rootfile + ".poli", region )
    dummy,dummy,pa = dump1( rootfile + ".pa", region )
    panew = paflip(pa) 
    dummy,dummy,paerr = dump1( rootfile + ".paerr", region )
    dummy,dummy,V = dump1( rootfile + ".V.cm", region )
    fout = open( outfile, "wb" )     # binary because I'll be pickling to it
    pickle.dump( [rootfile,region,chan,vlsr,I,poli,panew,paerr,V], fout )
    fout.close() 

# helper file for spectralPlot
def plotPolPanel( ax, label, vlsr, I, poli, V, ymin, ymax ) :
    ax.plot( vlsr, I, color='black', linewidth=0.5 )
    ax.plot( vlsr, poli, color='red', linewidth=0.5 )
    ax.plot( vlsr, V, color='blue', linewidth=0.5 )
    ax.set_xlim( -14., 23. )
    ax.tick_params( labelsize=10 )
    if ymin:
      ax.set_ylim( ymin, ymax )
    ax.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
    ax.annotate( "%s" % label, xy=(.95,.85), xycoords='axes fraction', \
        horizontalalignment="right", size=10 )

def plotPApanel( ax, vlsr10, pa10, pa10err, vlsr21, pa21, pa21err ) :
    ax.errorbar(vlsr10,pa10,yerr=pa10err,fmt='o',color='red') 
    ax.errorbar(vlsr21,pa21,yerr=pa21err,fmt='D',color='blue') 
    ax.set_xlim( -14., 22. )
    ax.set_ylim( 50., 130. )
    ax.grid( True, linewidth=0.02, color="0.02" )   # color=0.1 is a light gray
     
# specialized plotting routine for polarization plots
def polPlot( posList ) :
    pyplot.ioff()
    pp = PdfPages("polPlot.pdf")
    fig = pyplot.figure( figsize=(8,11) )
    xa = .1
    ya = .82

    for pos in posList :
      fin = open( pos["pos21"], "rb" ) 
      [rootfile,region,chan,vlsr21,I21,poli21,pa21,pa21err,V21] = pickle.load( fin )
      fin.close()
      ax = fig.add_axes( [xa,ya,.38,.15] )
      plotPolPanel( ax, "2-1", vlsr21, I21, poli21, V21, pos["ymin21"], pos["ymax21"]  )
      ax.set_xticklabels( [] )
      ax.annotate( "%s" % pos["label"], xy=(.5,.82), xycoords='axes fraction', \
        horizontalalignment="center", size=16 )

      fin = open( pos["pos10"], "rb" ) 
      [rootfile,region,chan,vlsr10,I10,poli10,pa10,pa10err,V10] = pickle.load( fin )
      fin.close()
      ya = ya - .15
      ax = fig.add_axes( [xa,ya,.38,.15] )
      plotPolPanel( ax, "1-0", vlsr10, I10, poli10, V10, pos["ymin10"], pos["ymax10"]  )

      #ya = ya - .1
      #ax = fig.add_axes( [xa,ya,.38,.1] )
      #plotPApanel( ax, vlsr10, pa10, pa10err, vlsr21, pa21, pa21err )
      ya = ya - .2
      if ya < .3 :
        xa = xa + .46
        ya = .82
    
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

def polPlot2( posList, j21=False ) :
    pyplot.ioff()
    pp = PdfPages("polPlot.pdf")
    fig = pyplot.figure( figsize=(8,11) )
    xa = .1
    ya = .78

    for pos in posList :
      if j21 :
        fin = open( pos["pos21"], "rb" ) 
        [rootfile,region,chan,vlsr21,I21,poli21,pa21,pa21err,V21] = pickle.load( fin )
        fin.close()
        ax = fig.add_axes( [xa,ya,.78,.2] )
        plotPolPanel( ax, "2-1", vlsr21, I21, poli21, V21, pos["ymin21"], pos["ymax21"]  )
        ax.set_xticklabels( [] )
        ax.annotate( "%s" % pos["label"], xy=(.5,.82), xycoords='axes fraction', \
          horizontalalignment="center", size=16 )

      else:
        fin = open( pos["pos10"], "rb" ) 
        [rootfile,region,chan,vlsr10,I10,poli10,pa10,pa10err,V10] = pickle.load( fin )
        fin.close()
        ax = fig.add_axes( [xa,ya,.78,.2] )
        plotPolPanel( ax, "1-0", vlsr10, I10, poli10, V10, pos["ymin10"], pos["ymax10"]  )
        ax.annotate( "%s" % pos["label"], xy=(.5,.82), xycoords='axes fraction', \
          horizontalalignment="center", size=16 )

      #ya = ya - .1
      #ax = fig.add_axes( [xa,ya,.38,.1] )
      #plotPApanel( ax, vlsr10, pa10, pa10err, vlsr21, pa21, pa21err )
      ya = ya - .22
      if ya < .05 :
        xa = xa + .46
        ya = .82
    
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

    
# compute percentiles of P/I vs I, V/I vs I, V/I vs P (NOT as function of velocity)
#def polStats2( IQUVmapList, pcutoff, icutoff, region='arcsec,box(1.2)', chanrange=[40,105] ):
def polStats2( IQUVmapList, pcutoff, icutoff, region='arcsec,box(1.2)', chanrange=[40,105] ):
    
  # write header info into each of the 3 output files
    fout = open("polStats2.dat","w")
    fout.write("# created by ori3.polStats2\n")
    fout.write("# input files: %s\n" % IQUVmapList[0] )
    fout.write("#              %s\n" % IQUVmapList[1] )
    fout.write("#              %s\n" % IQUVmapList[2] )
    fout.write("#              %s\n" % IQUVmapList[3] )
    fout.write("# region: %s\n" % region )
    fout.write("# chanrange: [%d,%d]\n" % (chanrange[0],chanrange[1]))
    fout.write("# icutoff: %.5f\n" % icutoff )
    fout.write("#\n#  Imin  Imax    25_pct 50_pct 75_pct 95_pct\n")

  # create giant array of I,Q,U,V for every pixel

    ifl = []
    fp = []
    fv = []
    vp = []
    z = [ [],[],[],[] ]    

    for ichan in range( chanrange[0],chanrange[1]+1) :
      print "processing chan %d" % ichan
      for nmap in range( 0,4 ):
        p = subprocess.Popen( ( shlex.split("imtab in=%s region=%s(%d) log=imtablog format=(3F12.5)" % 
           (IQUVmapList[nmap],region,ichan) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
        result = p.communicate()[0]
        x,y,flx = numpy.loadtxt("imtablog", unpack=True )
        z[nmap].append(flx)

  # now concatenate lists of numpy arrays into single arrays
    iarr = numpy.concatenate(z[0])
    qarr = numpy.concatenate(z[1])
    uarr = numpy.concatenate(z[2])
    varr = numpy.concatenate(z[3])

    deltaI = .05
    PI = []
    ABSV = []

    for Imin in numpy.arange(icutoff,2.,deltaI):
      print "Imin = %.5f" % Imin
      fracp = []
      fracv = []
      npts = 0
   
    # write percentiles of poli/I and abs(V)/I vs I
      for i,q,u,v in zip( iarr,qarr,uarr,varr ) :
        if (i > Imin) and (i <= (Imin+deltaI)) :
          npts = npts + 1
          poli = math.sqrt( q*q + u*u)
          fracp.append( poli/i )
          fracv.append( abs(v)/i )

        # fill out PI and ABSV arrays as we go
          PI.append(poli)
          ABSV.append( abs(v) )

      if npts > 50 :
        a = numpy.percentile( fracp, [25,50,75,95] )
        fout.write("%7.4f %7.4f %7.4f  %6d  %7.4f %7.4f %7.4f %7.4f" % \
          (Imin, Imin+deltaI, Imin+deltaI/2., npts, a[0],a[1],a[2],a[3]) )
        a = numpy.percentile( fracv, [25,50,75,95] )
        fout.write("    %7.4f %7.4f %7.4f %7.4f\n" % \
          (a[0],a[1],a[2],a[3]) )
    fout.close()

  # now use the accumulated PI and ABSV arrays to write percentiles of abs(V)/P vs P
  # note that all these points fulfil the requirement that i > Imin

    fout = open("polStats2.voverp","w")
    fout.write("# created by ori3.polStats2\n")
    fout.write("# input files: %s\n" % IQUVmapList[0] )
    fout.write("#              %s\n" % IQUVmapList[1] )
    fout.write("#              %s\n" % IQUVmapList[2] )
    fout.write("#              %s\n" % IQUVmapList[3] )
    fout.write("# region: %s\n" % region )
    fout.write("# chanrange: [%d,%d]\n" % (chanrange[0],chanrange[1]))
    fout.write("# icutoff: %.5f\n" % icutoff )
    fout.write("#\n#  Pmin  Pmax  Pmid  25_pct 50_pct 75_pct 95_pct\n")

    deltaP = .05
    for Pmin in numpy.arange(pcutoff,2.,deltaP):
      print "Pmin = %.5f" % Pmin
      voverp = []
      npts = 0
      for poli,absv in zip( PI, ABSV ) :
        if (poli > Pmin) and (poli <= (Pmin+deltaP)) :  
          voverp.append( absv/poli )
          npts = npts + 1
      if npts > 50 :
        a = numpy.percentile( voverp, [25,50,75,95] )
        fout.write("%7.4f %7.4f %7.4f  %6d  %7.4f %7.4f %7.4f %7.4f\n" % \
          (Pmin, Pmin+deltaP, Pmin+deltaP/2., npts, a[0],a[1],a[2],a[3]) )
    fout.close()

  # plot fracp, fracv vs I
  #  fig,ax = pyplot.subplots( nrows=2 )
  #  ax[0].scatter( ifl, fp, alpha=0.1 )
  #  pyplot.title("P/I vs I")
  #  ax[1].scatter( ifl, fv, alpha=0.1 )
  #  pyplot.title("V/I vs I")
  #  pyplot.show()
    

# retrieve restfreq and channel velocities from miriad map
def getVelInfo( miriadMap ) :
    p= subprocess.Popen( ( shlex.split('imlist in=%s' % miriadMap) ), \
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
    return restfreq,v1,dv

# panelPlot plots panels of pol data for a spectral line data cube
# Stokes V is in color; Stokes I is contours
# region is spatial region -- e.g., "arcsec,box(1)"
# chanRange gives images to plot -- e.g., [16,50]
# Vrange sets scale for Stokes V in Jy/beam -- e.g., [-.04,.04]
# contourList is for Stokes I in Jy/beam -- e.g., [.05,.1,.2,.4]

def panelPlot( IQUVmapList, region, chanRange, Vrange, contourList, outfile="panels.pdf" ) :
    delta = 4
    Pcutoff = .003
    Icutoff = .01
    vecScale = .3

    pp = PdfPages( outfile )
    ncols = 3
    nrows = 5 
    fig = pyplot.figure( figsize=(8,11) )
 
  # use ImageGrid to make neat array of boxes
    ax = ImageGrid( fig, 111, nrows_ncols=(nrows,ncols), share_all=True, axes_pad=0.1, \
      cbar_location="top", cbar_mode="single", cbar_pad=.3, cbar_size="3%" )

  # retrieve velocity axis
    restfreq,v1,dv = getVelInfo( IQUVmapList[0] )
  
  # read I,Q,U,V for each channel, plot 1 subpanel
    ichan = chanRange[0]  
    for ipanel in range(0, nrows*ncols ) :
      print "processing channel %d" % ichan 
      z = [ [],[],[],[] ]    
      for nmap in range( 0,4 ):
        p = subprocess.Popen( ( shlex.split("imtab in=%s region=%s(%d) log=imtablog format=(3F12.5)" % 
           (IQUVmapList[nmap],region,ichan) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
        result = p.communicate()[0]
        x,y,flx = numpy.loadtxt("imtablog", unpack=True )
        z[nmap].append(flx)
      nx = len( numpy.unique( x ) )
      ny = len( numpy.unique( y ) )
      print "nx = %d, ny = %d" % (nx,ny)
 
    # convert to 2d arrays
      I = numpy.reshape(z[0], (ny,nx))
      Q = numpy.reshape(z[1], (ny,nx))
      U = numpy.reshape(z[2], (ny,nx))
      V = numpy.reshape(z[3], (ny,nx))
      X = numpy.reshape(x, (ny,nx))
      Y = numpy.reshape(y, (ny,nx))

    # pixel plot of Stokes V
      im = ax[ipanel].imshow(V, origin='lower', aspect='equal', \
        extent=[x[0],x[-1],y[0],y[-1]], vmin=Vrange[0], vmax=Vrange[1], cmap='RdBu_r' )

    # set limits so that RA increases to the left
      ax[ipanel].set_xlim( x[0],x[-1] )
      ax[ipanel].set_ylim( y[0],y[-1] )

    # contour plot of Stokes I
      #ax[ipanel].contour( x.reshape(ny,nx), y.reshape(ny,nx), numpy.array(z[0]).reshape(ny,nx), \
      ax[ipanel].contour( X, Y, I, \
        colors='black', levels=contourList, linewidths=.4 )

    # compute poli and pa
      poli = numpy.sqrt( numpy.power(Q,2) + numpy.power(U,2) )
      pa = numpy.arctan2(Q,U)

    # plot line segment wherever poli > Pcutoff and I > Icutoff
    # use trick of creating one long line with None breaking up segments
    #  xlist = []
    #  ylist = []
    #  for n in range(0,nx,delta) :
    #    for m in range(0,ny,delta) :
    #      if (I[m,n] > Icutoff) and (poli[m,n] > Pcutoff) :
    #        lvec = poli[m,n]/I[m,n] * vecScale 
    #        xlist.append( X[m,n] - lvec/2. * math.sin(pa[m,n]) )
    #        xlist.append( X[m,n] + lvec/2. * math.sin(pa[m,n]) )
    #        xlist.append( None )
    #        ylist.append( Y[m,n] + lvec/2. * math.cos(pa[m,n]) )
    #        ylist.append( Y[m,n] - lvec/2. * math.cos(pa[m,n]) )
    #        ylist.append( None )
    #  if len(xlist) > 1 :
    #    ax[ipanel].plot( xlist, ylist, linestyle='-', color='white', linewidth=1)
    #    ax[ipanel].plot( xlist, ylist, linestyle='-', color='black', linewidth=0.5)

    # alternative method, coloring the vectors
    #  for n in range(0,nx,delta) :
    #    for m in range(0,ny,delta) :
    #      if (I[m,n] > Icutoff) and (poli[m,n] > Pcutoff) :
    #        xlist = [0,0]
    #        ylist = [0,0]
    #        lvec = poli[m,n]/I[m,n] * vecScale 
    #        xlist[0] = X[m,n] - lvec/2. * math.sin(pa[m,n]) 
    #        xlist[1] = X[m,n] + lvec/2. * math.sin(pa[m,n]) 
    #        ylist[0] = Y[m,n] + lvec/2. * math.cos(pa[m,n]) 
    #        ylist[1] = Y[m,n] - lvec/2. * math.cos(pa[m,n]) 
    #        color = pa[m,n]/math.pi + 0.5
    #        ax[ipanel].plot( xlist, ylist, linestyle='-', color='white', linewidth=1)
    #        ax[ipanel].plot( xlist, ylist, linestyle='-', color=cm.hsv(color), linewidth=0.5)

    # print velocity in upper right corner
      vlsr = v1 + dv*(ichan-1)
      ax[ipanel].text( .95, .95, "%5.1f" % vlsr, horizontalalignment="right", verticalalignment="top", \
         transform=ax[ipanel].transAxes )
      ax[ipanel].tick_params( axis='both', which='major', labelsize=8 )
      ichan = ichan + 1

  # finish up   
    ax.cbar_axes[0].colorbar(im)
    ax.cbar_axes[0].tick_params( axis='both', which='major', labelsize=8 )
    ax.cbar_axes[0].set_title("Stokes V (Jy/beam)\n", fontdict={'fontsize':8}) 
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()
