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

# vis1 is the leakage calibration file
vis1 = { "fileName"   : "3c273.pb",
         "selectStr"  : "source(3c273),time(05:00,23:59)",
         "optionStr"  : "circular,noxy,nopass",
         "flux"       : "4.04",
         "refant"     : 8, 
         "avgchan"    : 6,
         "lkname"     : "lkc" }

#vis2 = { "fileName"   : "3c279.pb",
#         "selectStr"  : "source(3c279)",
#         "optionStr"  : "circular,noxy,nopass",
#         "flux"       : "",
#         "refant"     : 8, 
#         "avgchan"    : 6,
#         "lkname"     : "dummy" }

vis2 = { "fileName"   : "SgrA.pb",
         "selectStr"  : "source(SGRA),uvrange(20,1000)",
         "optionStr"  : "circular,noxy,nopass",
         "flux"       : "",
         "refant"     : 8, 
         "avgchan"    : 6,
         "lkname"     : "dummy" }

     
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
    [src, I, rmsI, ncorrs] = paPlot.parseLine( lines[nlines-6] )
    [dummy, Q, rmsQ, ncorrs] = paPlot.parseLine( lines[nlines-5] )
    [dummy, U, rmsU, ncorrs] = paPlot.parseLine( lines[nlines-4] )
    [dummy, V, rmsV, ncorrs] = paPlot.parseLine( lines[nlines-3] )
  print "\n%7.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f  %6.4f %5.4f" % (I, rmsI, Q, rmsQ, U, rmsU, V, rmsV) 
  return [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ]


# Gpcopy leakage item from file1 to file2
def gpcopyLeak( file1, file2 ) :
  p= subprocess.Popen( shlex.split('gpcopy vis=%s out=%s options=nocal,nopass' \
     % (file1, file2 ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  print result


# Copy leakage item saved as file1 to data directory file2
def restoreLeak( file1, file2 ) :
  p= subprocess.Popen( shlex.split('cp -v %s %s/leakage' % (file1, file2 ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  print result


# Save leakage item from file1/leakage as file2
def saveLeak( file1, file2 ) :
  p= subprocess.Popen( shlex.split('cp -v %s/leakage %s' % (file1, file2 ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  print result

def makeLeakTables( ) :
  leakSolve.makeFtable( vis1 )
  for nwin in range (0,8) :
    for ch1 in range( 24*nwin+1, 24*nwin+24, 6 ) :
      ch2 = ch1 + vis2["avgchan"] - 1
      lineStrLeak = "chan,1,%d,%d" % ( ch1, vis2["avgchan"] )
      leakSolve.delHd( vis1 )
      leakSolve.addLeak( vis1, ch1 )
      saveLeak( vis1["fileName"], "leak.ch%03d-%03d" % (ch1,ch2) )


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


      #selectString = "source(0958+655),-ant(7)"
      #[ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ] = getStokes2( vis2["fileName"], selectString, lineStr )
      #fout = open("0958resultsB", "a")
      #fout.write("%8.3f %8.3f  %11.4e %11.4e  %11.4e %11.4e  %11.4e %11.4e  %11.4e %11.4e\n" \
      #   % (fstart,fstop,I,rmsI,Q,rmsQ,U,rmsU,V,rmsV) )
      #fout.close()

# Generate table of I,Q,U,V in file vis2["fileName"]
# Uses SAVED leakages in files of "leakName.chxxx-yyy"

def makeIQUtable( outfile, leakName ) :
  leakSolve.makeFtable( vis1 )
  leakSolve.makeFtable( vis2 )
  print " "
  print "LEAK: " + numpy.array_str( numpy.array(vis1["chfreq"])[0:5], precision=3, max_line_width=200 ) + " .... " \
     + numpy.array_str( numpy.array(vis1["chfreq"])[-5:], precision=3, max_line_width=200 )
  print "DATA: " + numpy.array_str( numpy.array(vis2["chfreq"])[0:5], precision=3, max_line_width=200 ) + " .... " \
     + numpy.array_str( numpy.array(vis2["chfreq"])[-5:], precision=3, max_line_width=200 )
  print " "

  fout = open( outfile, "w" )
  fout.write("# leaks from %s, data from %s\n" % ( vis1["fileName"], vis2["fileName"] ) )
  fout.write("#      f1     f2     I     rmsI     Q     rmsQ      U     rmsU     V    rmsV\n" )
  fout.close()

  for nwin in range (0,8) :
    for ch1 in range( 24*nwin+1, 24*nwin+24, 6 ) :

    # restore leakage from saved file
      ch2 = ch1 + vis2["avgchan"] - 1
      lineStrLeak = "chan,1,%d,%d" % ( ch1, vis2["avgchan"] )
      fstartLeak = vis1["chfreq"][ch1-1] - 0.5*vis1["chwidth"][ch1-1]
      fstopLeak = vis1["chfreq"][ch2-1] + 0.5*vis1["chwidth"][ch2-1]
      if leakName == None :
        leakSolve.delHd( vis1 )
      else :
        restoreLeak( "%s.ch%03d-%03d" % (leakName,ch1,ch2), vis2["fileName"] )

      lineStr = "chan,1,%d,%d" % ( ch1, vis2["avgchan"] )
      fstart = vis2["chfreq"][ch1-1] - 0.5*vis2["chwidth"][ch1-1]
      fstop = vis2["chfreq"][ch2-1] + 0.5*vis2["chwidth"][ch2-1]
      print "Leak:  %s   %8.3f %8.3f" % (lineStrLeak, fstartLeak, fstopLeak)
      print "Data:  %s   %8.3f %8.3f" % (lineStr, fstart, fstop)
      [ I, rmsI, Q, rmsQ, U, rmsU, V, rmsV ] = getStokes2( vis2["fileName"], vis2["selectStr"], lineStr )
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
  
  
