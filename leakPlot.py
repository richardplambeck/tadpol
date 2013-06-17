# leakPlot.py
#
#import subarrayCommands as SAC
import math
import time
#import device
import cmath
import numpy
import pylab
import sys

# leakage tables:
#  - separate table for each ant; default name is "lkxx" or "leakxx"
#  - each line: "f1 f2 DRamp DRphs DLamp DLphs DRre DRim DLre DLim pctQ pctU lineString"

# fileList:
#  - each line lists one directory and optional comment, e.g.
#    "/fringe2/plambeck/pol/leakage/08aug2012/lk" "3C279, 10 MHz steps" 
#  - program looks for files of form "dir%d % ant"

# types: 1 = amp vs freq

# read leaktables into memory, plot leakages of selected antennas from freq f1-f2
def leakPlot( fileList, ant, f1=0., f2=0. ) :
  files = []
  notes = []
  fin = open( fileList, "r" )
  for line in fin :
    print line
    if (len(line) > 0) and (line[0] != "#") :
      a = line.split()
      files.append(a[0])
      if len(a) > 1 :
        notes.append(a[1])
      else :
        notes.append(" ")
      for n in range(ant,ant+1) :
        try :
          filename = a[0]+"%d" % n
          din = open( filename, "r" )
        except IOerror :
          print "... can't open file %s" % filename
        else :
          print "... reading data from file %s" % filename
          leaklines = []
          for line in din :
            leaklines.append[line]   
          din.close()
           

# --- make vector plot of phased antennas or all antennas --- #
def plotVec( antList=phasedAnts ) :
  color = ['r','g','b','m','c','r','g','b','m','c','c','m','c','m','c']
  pbCorrectedVis = numpy.load("pbCorrectedVis.npy")
  #print pbCorrectedVis
  scalarSum = 0.
  vecList = [0.+0.j]	
  for n in antList :
    for m in range (0,16) :
      ilast = len(vecList) - 1
      if ( numpy.isnan( numpy.real( pbCorrectedVis[m][n-1] ) ) or \
           numpy.isnan( numpy.imag( pbCorrectedVis[m][n-1] ) ) ) :
        vecList = numpy.append( vecList, vecList[ilast] )
      else : 
        vecList = numpy.append( vecList, vecList[ilast] + pbCorrectedVis[m][n] )
        scalarSum += numpy.abs( pbCorrectedVis[m][n] )
  print "\nvecList:"
  print numpy.array_str( vecList, precision=2, max_line_width=200 )
  istart = 0
  i = 0
  while (istart < len(vecList) ) :
    x = numpy.real( vecList[istart:(istart+16)] )
    y = numpy.imag( vecList[istart:(istart+16)] )
    pylab.plot( x, y, color=color[i] )
    i = i+1
    istart = istart+16
  pylab.axis( [-scalarSum,scalarSum,-scalarSum,scalarSum] )
  pylab.grid(True)
  pylab.axes().set_aspect('equal')
  pylab.draw()

def plotVec2( antList=[1,2,3,4,5,6] ) :
  color = ['b','g','r','c','m','y','b','g','r','c','m','y','b','g','r']
  pbCorrectedVis = numpy.load("lastVis.npy")
  tsys = numpy.load("lastTsys.npy")
  snrMatrix = pbcorrectedVis * tsysMatrix * KperJyMatrix
	# element-wise multiplication of elements
  vecList = numpy.zeros( [15,17], dtype=complex )
	# ... list of (x,y) values, as complex numbers, for each vector, beginning with (0,0)
  maxScalarSum = 0. 
  pylab.ion()
  pylab.clf()
  for n in range(0,15) :
    scalarSum = 0.
    for m in range (0,16) :
      if ( numpy.isnan( numpy.real( pbCorrectedVis[m][n] ) ) or \
           numpy.isnan( numpy.imag( pbCorrectedVis[m][n] ) ) ) :
        vecList[n][m+1] = vecList[n][m]    # if nan, repeat last value
      else : 
        vecList[n][m+1] = vecList[n][m] + pbCorrectedVis[m][n] 
        scalarSum += numpy.abs( pbCorrectedVis[m][n] )
    if (scalarSum > maxScalarSum) :
      maxScalarSum = scalarSum
    print "\nvecList for ant %d:" % (n+1)
    print numpy.array_str( vecList[n], precision=2, max_line_width=200 )
  pylab.axis( [-1.05*maxScalarSum,1.05*maxScalarSum,-1.05*maxScalarSum,1.05*maxScalarSum] )
  for n in range(0,15) :
    x = numpy.real( vecList[n] )
    y = numpy.imag( vecList[n] )
    color = 'm'
    if (n+1) in antList :
      color = 'r'
    pylab.plot( x, y, color=color, linestyle='solid', linewidth=4 )
    pylab.annotate(str(n+1), [x[15],y[15]] )
  pylab.grid(True)
  pylab.axes().set_aspect('equal')
  pylab.draw()

def endless() :
  while True :
    plotVec( antList=phasedAnts )
    time.sleep(4.)

