# sefd.py
#
# this script computes sefd's for win13L,R and win15L,R for each ant

import numpy
import pylab
import sys
import os

jyPerK = numpy.empty( 15, dtype=float )
jyPerK[0:6] = 65.
jyPerK[6:15] = 145.

def calcSEFD( tsysfile, outfile ) :
  fout = open( outfile, "a" )
  fout.write("\n%s\n\n" % tsysfile)
  fin = open( tsysfile, "r" )
  for line in fin :
    if line.startswith("#") :
       sefd = tsys * jyPerK * gains * gains
       #print sefd
       sefdsort = numpy.sort(sefd)
       antsort = numpy.argsort(sefd)
       antsort = antsort + 1
       print antsort
       for m in range(0,15):
         if antsort[m] != 10 : 
           fout.write( "%3d  %5d  %5.2f  %6d\n" % (antsort[m], tsys[antsort[m]-1], gains[antsort[m]-1], sefdsort[m]) )
       fout.write( "\n")
    else :
      a = line.split()
      list = []
      for m in range(1,16) :
        list.append(float(a[m]))
      if ":" in a[0] :
        tsys = numpy.array( list )
        #print tsys
      else :
        gains = numpy.array( list )
        #print gains
  fout.close()
 
def doall( ):
  outfile = "sefdSummary.dat"
  calcSEFD ( "win13L.tsys", outfile )
  calcSEFD ( "win15L.tsys", outfile )
  calcSEFD ( "win13R.tsys", outfile )
  calcSEFD ( "win15R.tsys", outfile )
