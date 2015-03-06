# sefd.py
#
# this script computes sefd's for win13L,R and win15L,R for each ant

import vlbiCal
import numpy
import pylab
import sys
import os

jyPerK = numpy.empty( 15, dtype=float )
jyPerK[0:6] = 65.
jyPerK[6:15] = 145.

# inputs:
#    matrix of gains[times][ants]
#    matrix of tsys[times][ants]
# outputs:
#    [antsort,sefdsort] - sorted lists

def calcSEFD( gains, tsys ) :
  sefd = jyPerK * tsys * gains * gains   # element by element multiply, not matrix multiply
  #sefd = jyPerK * tsys   # element by element multiply, not matrix multiply
     # sefd is a matrix with dimension times
  sefdavg = numpy.average(sefd,axis=0)
  sefdsort = numpy.sort(sefdavg)
  antsort = numpy.argsort(sefdavg)
  antsort = antsort + 1
  return [antsort,sefdsort]


def doall( ):
  band = [ "5.L", "6.L", "7.L", "8.L", "5.R", "6.R", "7.R", "8.R" ] 
  antsort = numpy.zeros( (8,15), dtype=int )
  sefdsort = numpy.zeros( (8,15), dtype=float )
  # outfile = "sefdSummary.dat"
  [t, gainComplex] = vlbiCal.getGains( "wide.av" )
  gain = numpy.abs(gainComplex)
  print gain, len(gain)
  for i in range(0,8) :
    tsys = vlbiCal.getVar15( "tsys%s" % band[i], t )
    [antsort[i],sefdsort[i]] = calcSEFD( gain, tsys ) 
  print "%7s %15s %15s %15s %15s %15s %15s %15s" % (band[0],band[1],band[2],band[3],band[4],band[5],band[6],band[7])
  for n in range(0,15) :
    print "%4d %8.0f    %4d %8.0f    %4d %8.0f   %4d %8.0f    %4d %8.0f    %4d %8.0f   %4d %8.0f    %4d %8.0f" % \
      ( antsort[0][n],sefdsort[0][n], \
       antsort[1][n],sefdsort[1][n], \
       antsort[2][n],sefdsort[2][n], \
       antsort[3][n],sefdsort[3][n], \
       antsort[4][n],sefdsort[4][n], \
       antsort[5][n],sefdsort[5][n], \
       antsort[6][n],sefdsort[6][n], \
       antsort[7][n],sefdsort[7][n] )


