# ooLeak.py
#
# ori.py
# collection of routines to process ALMA data


import math
import time
import cmath
import numpy
import pylab
import sys
import RM
import pickle
import random
import matplotlib
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle

# the spectrum object - read spectral data created by casa viewer

class spec:
 
  def __init__(self, specfile ) :
    """read spectrum from disk file, sort by frequency"""
    amp = []
    freq = []
    self.file = specfile        # name of file in quotes
    try :
      fin = open( self.file, "r" )
    except :
      print "... can't open file %s" % self.file
    else :
      print "... reading data from file %s" % self.file

      for line in fin :
        if line.startswith("#") :
          print line[0:-1]
        elif len(line) > 1 :
          a = line.split()
          freq.append( float(a[0]) )
          amp.append( float(a[1])  )  
      fin.close()
   
    self.s  = sorted(zip(freq,amp))
    print "%d points in the spectrum" % len(self.s)

  def panel(self, p, fstart=0., fstop=1000.) :
    x = []
    y = []
    print fstart, fstop
    for q in self.s :
      if (q[0] > fstart) and (q[0] < fstop) :
        print q
        x.append( q[0] )
        y.append( q[1] )
    p.plot( x, y, color="r", linewidth=0.5 )

  def yminmax(self) :
    ymax = -1000.
    ymin = 1000.
    for q in self.s :
      if q[1] > ymax : ymax = q[1]
      if q[1] < ymin : ymin = q[1]
    diff = ymax - ymin
    return [ ymin - .05*diff, ymax + .05*diff + .3 ]

  def hist(self, fmin=0., fmax=0.) :
    pyplot.ioff()
    pp = PdfPages("spectrum.pdf")
    npanels = 1
    if fmin == 0 :
      npanels = 12
      fmin = self.s[0][0]     # first freq in spectrum
    if fmax == 0 :
      fmax = self.s[-1][0]	  # last freq in spectrum
    print fmin, fmax
    ymin,ymax = self.yminmax()
    
    if (npanels == 12) :
      ch1 = range(0,15360,3840/3)
    np = 0 
    for n in range(0,npanels) :
      np = np + 1    
      p = pyplot.subplot(3,1,np)
      p.grid( True, linewidth=0.1, color="0.05" )   # color=0.1 is a light gray
      fmin = self.s[ ch1[n] ][ 0 ] 
      fmax = self.s[ ch1[n] + 1279 ][ 0 ] 
      x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
      p.xaxis.set_major_formatter( x_formatter )
      p.tick_params( axis='both', which='major', labelsize=6 )
      p.axis( [fmin, fmax, ymin, ymax] )
      self.panel( p, fstart=fmin, fstop=fmax )
      self.ann( p, "linelist", fmin, fmax, ymax )
      if np % 3 == 0 :
        pyplot.savefig( pp, format='pdf', bbox_inches='tight' )
        np = 0
    pp.close()

  #def ann(self,p) :
  #  p.annotate("SO2 300K", xy=(649.81,-0.1), xytext=(649.81,.4), horizontalalignment="center", rotation="vertical", size=4, \
  #     arrowprops=dict( arrowstyle="->", linewidth=0.2 ) )

  def yvalue( self, fGHz ) :
    '''find amplitude of spectrum at freq fGHz - use for annotating lines'''
    ftest = self.s[0][0]
    yvalue = self.s[0][1]
    for q in self.s :
      if (q[0] < fGHz) and (q[0] > ftest)  : 
        ftest = q[0]
        yvalue = q[1]
    return yvalue
    
  def ann( self, p, linelistFile, fmin, fmax, ymax ) :
    fin = open( linelistFile, "r" )
    for line in fin :
      if (len(line) > 1) and (not line.startswith("#")) :
        a = line.split()
        linefreq = float(a[0])
        species = a[3]
        if (linefreq >= fmin) and (linefreq <= fmax) :
          yval = self.yvalue( linefreq )
          p.annotate( species, xy=(linefreq,yval), xytext=(linefreq,ymax), horizontalalignment="center", rotation="vertical", size=4, \
             arrowprops=dict( arrowstyle="->", linewidth=0.2 ) )

       

