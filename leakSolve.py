# leakSolve.py 
# 27jun2012 - begin adapting from paPlot.py
# 05may2013 - deleted unused routines paCalc, parseLine, getStokes; as far as I can 
#   see, this code has nothing in common with paPlot.py

import math
import time
import numpy
import subprocess
import shlex
import string


# dictionary files control operation
# .. fileName = visibility file containing leakage data
# .. selectStr = select string for gpcal - source, time, etc
# .. optionStr = option string for gpcal - with or without qusolve
# .. flux = I or I,Q,U
# .. refant
# .. lkname = output file name - code adds antenna number to the end!

vis1 = { "fileName"   : "wide.pb",
         "selectStr"  : "source(3c279),-ant(7)",
         "optionStr"  : "circular,noxy,nopass,qusolve",
         "avgchan"      :2,
         "flux"       : "10.",
         "refant"     : 8, 
         "lkname"     : "lk.3c279b." }

vis2 = { "fileName"   : "wide.pb",
         "selectStr"  : "source(3c279),-ant(7)",
         "optionStr"  : "circular,noxy,nopass,qusolve",
         "avgchan"      :47,
         "flux"       : "10.",
         "refant"     : 8, 
         "lkname"     : "lk.3c279c." }

vis3 = { "fileName"   : "wide.av",
         "selectStr"  : "source(3c279),-ant(7)",
         "optionStr"  : "circular,noxy,nopass,qusolve",
         "avgchan"      :1,
         "flux"       : "10.",
         "refant"     : 8, 
         "lkname"     : "lk.3c279d." }

# --- delete existing leakage file ---
def delHd( fileName ) :
  p= subprocess.Popen( ( shlex.split('delhd in=%s/leakage' % ( fileName ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  
# --- generate vis dictionary items "chfreq", "chwidth" - list of fstart,fstop for every correlator channel
def makeFtable( vis ) :
  nstart = []
  nchan = []
  fstart = []
  fstep = []
  chfreq = []
  chwidth = []
  p= subprocess.Popen( ( shlex.split('uvlist options=spectra vis=%s' % vis["fileName"] ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  uvlistOutput = p.communicate()[0]
  lines = uvlistOutput.split("\n")
  for line in lines :
    a = line.split()
    if "starting channel" in line :
      for n in range( 3, len(a) ) : nstart.append( int(a[n]) ) 
    if "number of channels" in line :
      for n in range( 4, len(a) ) : nchan.append( int(a[n]) )
    if "starting frequency" in line :
      for n in range( 3, len(a) ) : fstart.append( float(a[n]) ) 
    if "frequency interval" in line :
      for n in range( 3, len(a) ) : fstep.append( float(a[n]) )
  for n in range(0, len(nstart) ) :
    if nchan[n] > 0 :
      for i in range( 0, nchan[n] ) :
         chfreq.append(fstart[n] + i*fstep[n])
         chwidth.append(fstep[n])
  vis["nstart"] = nstart
  vis["nchan"] = nchan
  vis["chfreq"] = chfreq
  vis["chwidth"] = chwidth

# --- run gpcal on particular channel range, return src Q, U, and C1-C15 Dx,Dy complex leakages
def runGpcal( vis, lineString='', interval=5 ) :
  DR = numpy.zeros( 15, dtype=complex)
  DL = numpy.zeros( 15, dtype=complex)
  percentQ = 0.
  percentU = 0.
  p= subprocess.Popen( ( shlex.split('gpcal vis=%s select=%s line=%s flux=%s refant=%d interval=%0.1f options=%s' % \
     (vis["fileName"], vis["selectStr"], lineString, vis["flux"], vis["refant"], interval, vis["optionStr"] ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  gpcalOutput = p.communicate()[0]
  lines = gpcalOutput.split("\n")
  for line in lines :
    print line
    a = line.split()
    if "Percent Q:" in line :
      percentQ = float(a[2])
    if "Percent U:" in line :
      percentU = float(a[2])
    if "Dx,Dy" in line :
      ant = int(line[4:6])
      if ant < 16 :
        DR[ant-1] = float(line[16:22]) + 1.j * float(line[23:29])
        DL[ant-1] = float(line[32:38]) + 1.j * float(line[39:45])
  return [percentQ,percentU,DR,DL]

def addLeak( vis, ch1 ) :
  lineString = "chan,1,%d,%d" % ( ch1, vis["avgchan"] )
  ch2 = ch1 + vis["avgchan"] - 1
  fstart = vis["chfreq"][ch1-1] - 0.5*vis["chwidth"][ch1-1]
  fstop = vis["chfreq"][ch2-1] + 0.5*vis["chwidth"][ch2-1]
  print "%s   %8.3f %8.3f" % (lineString, fstart, fstop)
  [percentQ,percentU,DR,DL] = runGpcal( vis, lineString=lineString )
  for nant in range(0,15) :
    filename = vis["lkname"] + "%d" % (nant+1)
    fout = open( filename, "a" )
    fout.write("%8.3f %8.3f %8.3f %6.1f %8.3f %6.1f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f    %s\n" % \
       ( fstart, fstop, numpy.abs(DR[nant]), numpy.angle(DR[nant],deg=True), \
       numpy.abs(DL[nant]), numpy.angle(DL[nant],deg=True), DR[nant].real, \
       DR[nant].imag, DL[nant].real, DL[nant].imag, percentQ, percentU, lineString) )
    fout.close()

def addLk( vis, ch1 ) :
  lineString = "chan,1,%d,%d" % ( ch1, vis["avgchan"] )
  ch2 = ch1 + vis["avgchan"] - 1
  fstart = vis["chfreq"][ch1-1] - 0.5*vis["chwidth"][ch1-1]
  fstop = vis["chfreq"][ch2-1] + 0.5*vis["chwidth"][ch2-1]
  print "%s   %8.3f %8.3f" % (lineString, fstart, fstop)
  [percentQ,percentU,DR,DL] = runGpcal( vis, lineString=lineString )
  filename = vis["LkFile"]
  fout = open( filename, "a" )
  for nant in range(0,15) :
    fout.write("%3d %8.3f %8.3f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f    %s\n" % \
       ( nant+1, fstart, fstop, DR[nant].real, DR[nant].imag, DL[nant].real, \
       DL[nant].imag, percentQ, percentU, lineString) )
  fout.close()

def makeTable( vis ) :
  for nant in range(0,15) :
    filename = vis["lkname"] + "%d" % (nant+1)
    fout = open( filename, "w" )    # wipes out previous file!
    fout.write("# input data: %s\n" % vis["fileName"] )
    fout.write("# selectStr : %s\n" % vis["selectStr"] )
    fout.write("# optionStr : %s\n" % vis["optionStr"] )
    fout.write("# flux      : %s\n" % vis["flux"] )
    fout.write("# refant    : %d\n" % vis["refant"] )
    fout.write("# avgchan   : %d\n" % vis["avgchan"] )
    fout.close()
  makeFtable( vis )	   # fill in nstart, nchan, chfreq, chwidth dictionary items
  for nstart,nchan in zip( vis["nstart"],vis["nchan"] ) :
    if nchan > 0 :
      for n1 in range( nstart, nstart+nchan-vis["avgchan"]+1, vis["avgchan"] ) :
        print n1, n1+vis["avgchan"] - 1
        addLeak( vis, n1 )
  
def makeLk( vis ) :
  filename = vis["LkFile"]
  fout = open( filename, "w" )    # wipes out previous file!
  fout.write("# input data: %s\n" % vis["fileName"] )
  fout.write("# legend    : %s\n" % vis["legend"] )
  fout.write("# selectStr : %s\n" % vis["selectStr"] )
  fout.write("# optionStr : %s\n" % vis["optionStr"] )
  fout.write("# flux      : %s\n" % vis["flux"] )
  fout.write("# refant    : %d\n" % vis["refant"] )
  fout.write("# avgchan   : %d\n" % vis["avgchan"] )
  fout.close()
  makeFtable( vis )	   # fill in nstart, nchan, chfreq, chwidth dictionary items
  for nstart,nchan in zip( vis["nstart"],vis["nchan"] ) :
    if nchan > 0 :
      for n1 in range( nstart, nstart+nchan-vis["avgchan"]+1, vis["avgchan"] ) :
        print n1, n1+vis["avgchan"] - 1
        addLk( vis, n1 )

# copies Lk specified by lineStr from LkFile to outfile
def restoreLk( LkFile, lineStr, outfile ) :
  delHd( outfile )
  failed = True
  try:
    fin = open( LkFile, "r" )
  except :
    print "couldn't locate file %s" % LkFile
    return
  lk = numpy.zeros( 60 )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      if len(a) > 9 :
        if a[9] == lineStr :
          failed = False
          fstart = float(a[1])
          fstop = float(a[2])
          nant = int(a[0])
          lk[4*nant-4] = float(a[3])
          lk[4*nant-3] = float(a[4])
          lk[4*nant-2] = float(a[5])
          lk[4*nant-1] = float(a[6])
  if failed :
    print "FAILED: could not find %s in file %s" % (lineStr,LkFile)
  else :
    leakStr = ""
    for n in range(0,60) :
      if n == 0 :
        leakStr = "%.3f" % lk[0]
      else :
        leakStr = leakStr + ",%.3f" % (lk[n])
    print "leakStr=%s" % leakStr 
    p= subprocess.Popen( ( shlex.split('/o/plambeck/1mmDualPol/Calibration/writeleak vis=%s leak=%s' % \
     ( outfile, leakStr ) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
    result = p.communicate()[0]
    print result 
 
  
# ======================== start semi-obsolete wip routiones ========================= #

def makeWip( infile, outfile, ncolumn, hist=False, dots=False ) :
  prevx = 0.
  ncol = ncolumn - 1
  fin = open( infile, "r" )
  fout = open( outfile, "a" )
  for line in fin:
    a = line.split()
    if a[ncol] != "nan" :
      f1 = float(a[0])
      f2 = float(a[1])
      fmean = 0.5*(f1+f2)
      fwidth = abs(f2-f1)
      if not hist :          # smooth curve
        if abs(prevx - fmean) < (1.1*fwidth) :
          fout.write("draw %8.3f %s\n" % (fmean,a[ncol]) )
        else :
          fout.write("move %8.3f %s\n" % (fmean,a[ncol]) )
        if dots :
          fout.write("dot\n")  
        prevx = fmean
      else :                 # histogram
        if (abs(prevx-float(a[0])) < 0.1*fwidth)  and connect :
          fout.write("draw %s %s\n" % (a[0],a[ncol]) )
        else :
          fout.write("move %s %s\n" % (a[0],a[ncol]) )
        fout.write("draw %s %s\n" % (a[1],a[ncol]) ) 
      prevx = float(a[1])
  fin.close()
  fout.close()

def addLine( outfile, linetext ) :
  fout = open( outfile, "a")
  fout.write("%s\n" % linetext)
  fout.close()
  
def save() :
  for nant in range(1,16) :
    infile = "leak%d" % nant
    outfile = "DxAmp%d.wip" % nant
    makeWip( infile, outfile, 3)
    outfile = "DyAmp%d.wip" % nant
    makeWip( infile, outfile, 5)

    infile = "/fringe2/plambeck/pol/leakage/24jun2012/leak%d" % nant
    outfile = "DxAmp%d.wip" % nant
    addLine( outfile, "color 2")
    makeWip( infile, outfile, 3)
    outfile = "DyAmp%d.wip" % nant
    addLine( outfile, "color 2")
    makeWip( infile, outfile, 5)

# ======================== start semi-obsolete wip routiones ========================= #
  
def doit() :
  makeTable( vis3 ) 
    

  
  
