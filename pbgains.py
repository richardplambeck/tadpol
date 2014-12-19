# pbgains.py
#
# computes avg passband from gpplt log file
# measuring notch at 5 GHz - 11dec2014

import math
import time
import numpy
import subprocess
import shlex

antList = [2,3,4,5,8,12,13,14]
 
# read log file created by "gpplt vis=xx options=bandpass log=yy"

def read_gpplt_log (infile, avgfile ) :
  fin = open( infile, "r" )
  fout = open( avgfile, "w" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      if len(a) == 7 :
        freq = float( a[0] )
        g = []
        for n in range(1,7) :
          g.append( a[n] )
      elif len(a) == 6 :
        for n in range(0,6) :
          g.append( a[n] )
      elif len(a) == 5 :
        for n in range(0,5) :
          g.append( a[n] )
        gain = numpy.array( g, dtype=float )
        wgt = numpy.ones( 23 )
        for n in range(0,len(g)) :
          if gain[n] < 0.001 :
            wgt[n] = 0.
        print gain
        print wgt
        gavg = numpy.average( gain, weights=wgt )
        print freq, gavg
        fout.write("%10.5f  %.5f\n" % (freq-90., gavg) )  
  fin.close()
  fout.close()
