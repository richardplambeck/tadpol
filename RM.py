# RM.py
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

     
# Extract Stokes (I,Q,U, or V) amplitude and rms from one line of uvflux output 
def parseLine( line ) :
  print line
  a = line.split()
  if (len(a) == 9) :
    src = a[0]
    stokes = float(a[3])
    ncorrs = int(a[8])
    rms = float(a[5])/math.sqrt(ncorrs)
  else :
    src = ' '
    stokes = float(a[2])
    ncorrs = int(a[7])
    rms = float(a[4])/math.sqrt(ncorrs)
  return [ src, stokes, rms, ncorrs ]

# Run uvflux on selected data; return I, Q, U, V, and rms errors
def getStokes2( infile, selectString, lineString ) :
  p= subprocess.Popen( ( shlex.split('uvflux vis=%s select=%s line=%s stokes=I,Q,U,V' \
     % (infile, selectString, lineString) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  lines = result.split("\n")
  I = rmsI = pa = sigmaPa = poli = sigmaPoli = V = rmsV = 0.	# in case we get error
  Q = rmsQ = U = rmsU = V = rmsV = 0.
  nlines = len(lines)
  if nlines >= 10 :
    [src, I, rmsI, ncorrs] = parseLine( lines[nlines-6] )
    [dummy, Q, rmsQ, ncorrs] = parseLine( lines[nlines-5] )
    [dummy, U, rmsU, ncorrs] = parseLine( lines[nlines-4] )
    [dummy, V, rmsV, ncorrs] = parseLine( lines[nlines-3] )
  print "\n%7.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f" % (I, rmsI, Q, rmsQ, U, rmsU, V, rmsV) 
  return [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ]

# makes list of lineStrings in LkFile
def makeStrList( LkFile ) :
  fin = open( LkFile, "r" )
  strList = []
  for line in fin :
    if not line.startswith("#") and len(line) > 1 :
      a = line.split()
      if a[0] == "C01" :
        strList.append( a[9] )
  fin.close()
  return strList


# Generate table of I,Q,U,V in file vis["fileName"]
# Uses SAVED leakages in files of "leakName.chxxx-yyy"

def makeIQUtable( vis, outfile, LkFile ) :
  strList = makeStrList( LkFile )
  leakSolve.makeFtable( vis )

  #leakSolve.makeFtable( vis2 )
  #print " "
  #print "LEAK: " + numpy.array_str( numpy.array(vis1["chfreq"])[0:5], precision=3, max_line_width=200 ) + " .... " \
  #   + numpy.array_str( numpy.array(vis1["chfreq"])[-5:], precision=3, max_line_width=200 )
  #print "DATA: " + numpy.array_str( numpy.array(vis2["chfreq"])[0:5], precision=3, max_line_width=200 ) + " .... " \
  #   + numpy.array_str( numpy.array(vis2["chfreq"])[-5:], precision=3, max_line_width=200 )
  #print " "

  fout = open( outfile, "w" )
  fout.write("# leaks from %s, data from %s\n" % ( LkFile, vis["fileName"] ) )
  fout.write("#      f1     f2     I     rmsI     Q     rmsQ      U     rmsU     V    rmsV\n" )
  fout.close()

  for lineStr in strList :
    [fstartLeak,fstopLeak] = leakSolve.restoreLk( LkFile, lineStr, vis["fileName"] )
    a = lineStr.split(",")   # "chan,1,49,8"
    ch1 = int(a[2])
    ch2 = ch1 + int(a[3]) - 1
    fstart = vis["chfreq"][ch1-1] - 0.5*vis["chwidth"][ch1-1]
    fstop = vis["chfreq"][ch2-1] + 0.5*vis["chwidth"][ch2-1]
    print "  Leak %7.3f - %7.3f" % (fstartLeak, fstopLeak)
    print "  Data %7.3f - %7.3f" % (fstart, fstop)
    [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ] = getStokes2( vis["fileName"], vis["selectStr"], lineStr )
    fout = open( outfile, "a" )
    fout.write("%8.3f %8.3f   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e\n" \
         % (fstart,fstop,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV) )
    fout.close()

  
def RMsearch( infile, outfile ) :
  RMmin = -5.e9
  RMmax = 5.e9
  RMstep = 5.e5

# read in the data
  plist = []
  flist = []
  fin = open( infile, "r" )
  for line in fin :
    if not line.startswith("#") :
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
  
# ==== deprecated below this line ==== #  

# Make table of I,Q,U,V in file vis2["fileName"], using leakages from file vis1["fileName"]
def makeTable( outfile ) :
  leakSolve.makeFtable( vis1 )
  leakSolve.makeFtable( vis2 )
  print " =============== LEAK ============= "
  print numpy.array_str( numpy.array(vis1["chfreq"]), precision=3, max_line_width=200 )
  print " =============== DATA ============= "
  print numpy.array_str( numpy.array(vis2["chfreq"]), precision=3, max_line_width=200 )

  fout = open("outfile", "w")
  fout.write("# leaks from %s, data from %s\n" % ( vis1["fileName"], vis2["fileName"] ) )
  fout.write("#      f1     f2     I     rmsI     Q     rmsQ      U     rmsU     V    rmsV\n" )
  fout.close()

  for nwin in range (0,8) :
    for ch1 in range( 47*nwin+3,47*nwin+40,7) :

    # Solve for leakage
      ch2 = ch1 + vis2["avgchan"] - 1
      lineStrLeak = "chan,1,%d,%d" % ( ch1, vis2["avgchan"] )
      fstartLeak = vis1["chfreq"][ch1-1] - 0.5*vis1["chwidth"][ch1-1]
      fstopLeak = vis1["chfreq"][ch2-1] + 0.5*vis1["chwidth"][ch2-1]
      leakSolve.delHd( vis1 )
      leakSolve.addLeak( vis1, ch1 )
      gpcopyLeak( vis1["fileName"], vis2["fileName"] )
      saveLeak( vis1["fileName"], "leakC.ch%03d-%03d" % (ch1,ch2) )

      selectString = "source(M81),-ant(7)"
    # this is because Doppler tracking screwed everything up
      if (nwin < 4) :
        ch1 = ch1+1
        ch2 = ch2+1
      else :
        ch1 = ch1-1
        ch2 = ch2-1
      lineStr = "chan,1,%d,%d" % ( ch1, vis2["avgchan"] )
      fstart = vis2["chfreq"][ch1-1] - 0.5*vis2["chwidth"][ch1-1]
      fstop = vis2["chfreq"][ch2-1] + 0.5*vis2["chwidth"][ch2-1]
      print "Leak:  %s   %8.3f %8.3f" % (lineStrLeak, fstartLeak, fstopLeak)
      print "Data:  %s   %8.3f %8.3f" % (lineStr, fstart, fstop)
      [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ] = getStokes2( vis2["fileName"], selectString, lineStr )
      fout = open("M81resultsC", "a")
      fout.write("%8.3f %8.3f   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e   %10.3e %10.3e\n" \
         % (fstart,fstop,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV) )
      fout.close()

