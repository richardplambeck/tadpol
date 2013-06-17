# M81.py
#
# bits from paPlot.py, etc

import math
import time
import numpy
import subprocess
import shlex
import string
import leakSolve
import paPlot

# vis1 is the leakage calibration file
vis1 = { "fileName"   : "/fringe2/plambeck/pol/M81/18mar2013/wide.pb",
         "selectStr"  : "source(3c279),-ant(7)",
         "optionStr"  : "circular,noxy,nopass",
         "flux"       : "9.506,.452,.536",
         "refant"     : 8, 
         "avgchan"    : 5,
         "lkname"     : "lk." }

vis2 = { "fileName"   : "wide.cal",
         "selectStr"  : "source(M81),-ant(7)",
         "optionStr"  : "circular,noxy,nopass",
         "flux"       : "",
         "refant"     : 8, 
         "avgchan"    : 5,
         "lkname"     : "dummy" }

     
# run uvflux on selected data, return I, POLI, PA, and error estimates #
def getStokes2( infile, selectString, lineString ) :
  p= subprocess.Popen( ( shlex.split('uvflux vis=%s select=%s line=%s stokes=I,Q,U,V' \
     % (infile, selectString, lineString) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  lines = result.split("\n")
  I = rmsI = pa = sigmaPa = poli = sigmaPoli = V = rmsV = 0.	# in case we get error
  nlines = len(lines)
  if nlines >= 10 :
    [src, I, rmsI, ncorrs] = paPlot.parseLine( lines[nlines-6] )
    [dummy, Q, rmsQ, ncorrs] = paPlot.parseLine( lines[nlines-5] )
    [dummy, U, rmsU, ncorrs] = paPlot.parseLine( lines[nlines-4] )
    [dummy, V, rmsV, ncorrs] = paPlot.parseLine( lines[nlines-3] )
  print "%7.3f %5.3f  %6.3f %5.3f  %6.3f %5.3f  %6.3f %5.3f" % (I, rmsI, Q, rmsQ, U, rmsU, V, rmsV) 
  return [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ]

def runGpcopy( file1, file2 ) :
  p= subprocess.Popen( shlex.split('gpcopy vis=%s out=%s options=nocal,nopass' \
     % (file1, file2 ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  print result

def makeTable( ) :
  leakSolve.makeFtable( vis1 )
  leakSolve.makeFtable( vis2 )
  for f1,f2 in zip( vis1["chfreq"],vis2["chfreq"] ) :
    print f1-f2
  fout2 = open("0958results", "a")
  selectString = "source(M81),-ant(7)"
  for nwin in range (0,8) :
    for ch1 in range( 47*nwin+2,47*nwin+47,5) :
      ch2 = ch1 + vis2["avgchan"] - 1
      lineStr = "chan,1,%d,%d" % ( ch1, vis2["avgchan"] )
      fstart = vis2["chfreq"][ch1-1] - 0.5*vis2["chwidth"][ch1-1]
      fstop = vis2["chfreq"][ch2-1] + 0.5*vis2["chwidth"][ch2-1]
      print "%s   %8.3f %8.3f" % (lineStr, fstart, fstop)
      leakSolve.delHd( vis1 )
      leakSolve.addLeak( vis1, ch1 )
      runGpcopy( vis1["fileName"], vis2["fileName"] )

      selectString = "source(M81),-ant(7)"
      [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ] = getStokes2( vis2["fileName"], selectString, lineStr )
      fout = open("M81results", "a")
      fout.write("%8.3f %8.3f  %11.4e %11.4e  %11.4e %11.4e  %11.4e %11.4e  %11.4e %11.4e\n" \
         % (fstart,fstop,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV) )
      fout.close()


      selectString = "source(0958+655),-ant(7)"
      [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ] = getStokes2( vis2["fileName"], selectString, lineStr )
      fout = open("0958results", "a")
      fout.write("%8.3f %8.3f  %11.4e %11.4e  %11.4e %11.4e  %11.4e %11.4e  %11.4e %11.4e\n" \
         % (fstart,fstop,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV) )
      fout.close()

  
def RMsearch( infile, outfile ) :
  RMmin = -2.e9
  RMmax = 2.e9
  RMstep = 1.e6

# read in the data
  plist = []
  flist = []
  fin = open( infile, "r" )
  for line in fin :
    a = line.split()
    flist.append( (float(a[0]) + float(a[1]))/2. )
    Q = float(a[4])
    U = float(a[6])
    plist.append( Q + 1j*U )
  fin.close()  
  fGHz = numpy.array(flist)
  lambdasq = (.2998*.2998)/(fGHz*fGHz)
  pol = numpy.array(plist)
  norm = 1./len(pol)

  fout = open( outfile, "w" ) 
  for RM in range( RMmin, RMmax+RMstep, RMstep ) :
    rot = numpy.exp( 1j * RM * lambdasq )
    #print numpy.angle(rot, deg=True)
    P = norm * numpy.abs(numpy.sum( pol * rot )) 
    print "RM = %.3e, P = %5.3f" % (RM,P)
    fout.write( "%.3e  %11.4e\n" % ( RM, P ) )
  fout.close()
  
  
