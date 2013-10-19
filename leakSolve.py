# leakSolve.py 
# 27jun2012 - begin adapting from paPlot.py
# 05may2013 - deleted unused routines paCalc, parseLine, getStokes; as far as I can 
#   see, this code has nothing in common with paPlot.py
# 16jun2013 - delete code making separate leak files for each antenna - put all leakages into one Lk file
# 18oct2013 - add avgchan=0 option to average all channels in each window (this will be the default)

import math
import time
import numpy
import subprocess
import shlex
import string
import sys

# this is a collection of routines (not objects!) that solve for or restore leakages;
#    leakages are written to or read from a disk file
# 
# note that each correlator section will have a leakage; it is not possible to
#    average across the whole LSB or USB, or across all bands
#
# begin with a fully calibrated data set; it may include many sources, but only
#    a single source should be specified in selectStr
#
# dictionary file controls operation - generally this will be created in calling routine
# vis1 = { "fileName"   : "wide.cal",                    (no default)
#          "selectStr"  : "source(3c279),-ant(7)",       (no default)
#          "LkName"     : "Lk.22mar2013",                (no default)
#          "avgchan"    : 5,                             (default = 0 -> avg across corr section)
#          "interval"   : 5,                            
#          "refant"     : 8,                           
#          "optionStr"  : "circular,noxy,nopass,qusolve",
#          "flux"       : "10.",
#          "legend"     : "vlbi/22mar2013/3c279"          
#
# leakSolve.makeLk( vis )                                       - generates Lkfile (erase previous)
# leakSolve.restoreLk( LkFile, lineStr, outfile )               - restores leakage solution to outfile 
# leakSolve.stripout( Lkfile, antList, fstart, fstop, outfile ) - phased leakage sum for vlbi

# --- delete existing leakage file ---
def delHd( fileName ) :
  p= subprocess.Popen( ( shlex.split('delhd in=%s/leakage' % ( fileName ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  
# --- generate vis dictionary items "chfreq", "chwidth" - list of fstart,fstop for every correlator channel
def makeFtable( vis ) :
  'fills in dictionary items nstart,nchan (per section) and chfreq,chwidth (per chan)'
  nstart = []
  nchan = []
  fstart = []
  fstep = []
  chfreq = []
  chwidth = []
  p= subprocess.Popen( ( shlex.split('uvlist options=spectra vis=%s' % vis["fileName"] ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  uvlistOutput = p.communicate()[0]
  #print uvlistOutput
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
  #print "vis[nstart] = ", vis["nstart"]
  #print "vis[nchan] = ", vis["nchan"]
  #print "vis[chfreq] = ", vis["chfreq"]
  #print "vis[chwidth] = ", vis["chwidth"]
  #sys.stdout.flush()

# --- run gpcal on specified channel range, return src Q, U, and Dx,Dy complex leakages
def runGpcal( vis, lineString='' ) :
  DR = numpy.zeros( 15, dtype=complex)
  DL = numpy.zeros( 15, dtype=complex)
  percentQ = 0.
  percentU = 0.
  p= subprocess.Popen( ( shlex.split('gpcal vis=%s select=%s line=%s flux=%s refant=%d interval=%0.1f options=%s' % \
     (vis["fileName"], vis["selectStr"], lineString, vis["flux"], vis["refant"], vis["interval"], vis["optionStr"] ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  gpcalOutput = p.communicate()[0]
  lines = gpcalOutput.split("\n")
  for line in lines :
    a = line.split()
    if "Percent Q:" in line :
      percentQ = float(a[2])
    if "Percent U:" in line :
      percentU = float(a[2])
    if "Dx,Dy" in line :
      ant = int(line[4:6])
      if ant < 16 :
        print line
        DR[ant-1] = float(line[16:22]) + 1.j * float(line[23:29])
        DL[ant-1] = float(line[32:38]) + 1.j * float(line[39:45])
  return [ percentQ, percentU, DR, DL ]

# solve for leakages for a particular channel range, append to Lkfile
def addLk( vis, ch1, avgchan ) :
  lineString = "chan,1,%d,%d" % ( ch1, avgchan )
  ch2 = ch1 + avgchan - 1
  fstart = vis["chfreq"][ch1-1] - 0.5*vis["chwidth"][ch1-1]
  fstop = vis["chfreq"][ch2-1] + 0.5*vis["chwidth"][ch2-1]
  print ""
  print "%s   %8.3f %8.3f" % (lineString, fstart, fstop)
  [ percentQ, percentU, DR, DL ] = runGpcal( vis, lineString=lineString )
  filename = vis["LkFile"]
  fout = open( filename, "a" )
  fout.write("#\n")
  for nant in range(0,15) :
    fout.write("C%02d %8.3f %8.3f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f    %s\n" % \
       ( nant+1, fstart, fstop, DR[nant].real, DR[nant].imag, DL[nant].real, \
       DL[nant].imag, percentQ, percentU, lineString) )
  fout.close()

# fill in any missing vis dictionary items with defaults
def fillDefaults( vis ) :
  if "fileName" in vis :
    print "... fileName : ", vis["fileName"]
  else :
    print "... FATAL: no fileName specified"
  if not "selectStr" in vis :
    vis[ "selectStr" ] = ""  
    print "... WARNING: selectStr is blank - using all sources"
  print "... selectStr : ", vis[ "selectStr" ]
  if not "flux" in vis :
    vis["flux"] = "1."
  if not "interval" in vis :
    vis[ "interval" ] = 5.
  if not "LkName" in vis :
    vis[ "LkName" ] = "Lk"
  if not "refant" in vis :
    vis[ "refant" ] = 8
  if not "legend" in vis :
    vis[ "legend" ] = ""
    
# creates Lk file using prescription given in vis dictionary
def makeLk( vis ) :
  fillDefaults( vis )
  filename = vis["LkFile"]
  fout = open( filename, "w" )    # wipes out previous file!
  fout.write("# input data : %s\n" % vis["fileName"] )
  fout.write("# legend     : %s\n" % vis["legend"] )
  fout.write("# selectStr  : %s\n" % vis["selectStr"] )
  fout.write("# optionStr  : %s\n" % vis["optionStr"] )
  fout.write("# flux       : %s\n" % vis["flux"] )
  fout.write("# refant     : %d\n" % vis["refant"] )
  fout.write("# avgchan    : %d\n" % vis["avgchan"] )
  fout.write("# interval   : %.2f\n" % vis["interval"] )
  fout.close()
  makeFtable( vis )	   # fill in nstart, nchan, chfreq, chwidth dictionary items
  for nstart,nchan in zip( vis["nstart"],vis["nchan"] ) :
    if nchan > 0 :
      avgchan = vis["avgchan"]
      if (avgchan > nchan) or (avgchan <= 0) :
        avgchan = nchan
      nextra = nchan - (nchan/avgchan)*avgchan
      for n1 in range( nstart + nextra/2, nstart+nchan-avgchan+1, avgchan ) :
        print n1, n1+avgchan - 1
        delHd( vis["fileName"] ) 
        addLk( vis, n1, avgchan )

# copies Lk specified by lineStr from LkFile to outfile
def restoreLk( LkFile, lineStr, outfile ) :
  print "   restoreLk: copy leakages for lineStr = %s from %s to %s" % (lineStr, LkFile, outfile)
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
          nant = int(a[0].lstrip("C"))    # e.g., "C01" -> 1
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
    #print "leakStr=%s" % leakStr 
    p= subprocess.Popen( ( shlex.split('/fringe2/plambeck/pol/tadpol/writeleak vis=%s leak=%s' % \
     ( outfile, leakStr ) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
    result = p.communicate()[0]
    if (len(result) > 1) : print result 
  return [fstart,fstop]

# make table of D_R and D_L for the phased sum of ants in antList (or for single ant)
# over the freq range f1-f2, and produce the vector avg over that freq interval
def stripout( Lkfile, antList, fstart, fstop, outfile ) :
  fin = open( Lkfile, "r" )
  fout = open( outfile, "a" )
  fout.write("# LkFile: %s\n" % Lkfile )
  fout.write("# antList: %s\n" % str(antList) )
  fout.write("# fstart,fstop: %.3f,%.3f\n" % (fstart,fstop) )
  navg = 0
  ngrandavg = 0
  DRgrandavg = 0. + 1j * 0.
  DLgrandavg = 0. + 1j * 0.
  for line in fin :
    if (len(line) > 2) and line.startswith("#") :
      print line[0:-1]
      fout.write( "%s" % line )
    elif (len(line) > 2) :
      a = line.split() 
      ant = int( a[0].strip("C") )
      if ant == 1 :      # new freq interval
        if (navg > 0) :
          DRavg = DRavg/navg
          DLavg = DLavg/navg
          fout.write(" %8.3f %8.3f  %8.3f %8.3f %6.3f   %8.3f %8.3f %6.3f    %s\n" % \
            (f1, f2, abs(DRavg), DRavg.real, DRavg.imag, abs(DLavg), DLavg.real, \
            DLavg.imag, a[9] ) )
          DRgrandavg = DRgrandavg + DRavg
          DLgrandavg = DLgrandavg + DLavg
          ngrandavg = ngrandavg + 1
          print ngrandavg, DRgrandavg, DLgrandavg
        f1 = float(a[1])
        f2 = float(a[2])
        navg = 0
        DRavg = 0. + 1j * 0.
        DLavg = 0. + 1j * 0.
      if (ant in antList) and ((( f1 >= fstart ) and (f1 <= fstop)) or ((f2 >= fstart) and (f2 <= fstop))) :
        DR = float(a[3]) + 1j * float(a[4])
        DL = float(a[5]) + 1j * float(a[6])
        if (abs(DR) == 0.) and (abs(DL) == 0.) :
          print "entry for antenna %d is zero" % ant
          fout.write("  --- no solution for antenna %d ---\n" % ant )
        else :
          navg = navg + 1
          DRavg = DRavg + DR
          DLavg = DLavg + DL
  # finish up
  if (navg > 0) :
    DRavg = DRavg/navg
    DLavg = DLavg/navg
    fout.write(" %8.3f %8.3f  %8.3f %8.3f %6.3f   %8.3f %8.3f %6.3f    %s\n" % \
      (f1, f2, abs(DRavg), DRavg.real, DRavg.imag, abs(DLavg), DLavg.real, \
      DLavg.imag, a[9] ) )
    DRgrandavg = DRgrandavg + DRavg
    DLgrandavg = DLgrandavg + DLavg
    print ngrandavg, DRgrandavg, DLgrandavg
    ngrandavg = ngrandavg+1
  if ngrandavg > 0 :
    DRavg = DRgrandavg/ngrandavg
    DLavg = DLgrandavg/ngrandavg
    fout.write("\n  band avg          %8.3f %8.3f %6.3f   %8.3f %8.3f %6.3f \n\n\n" % \
      (abs(DRavg), DRavg.real, DRavg.imag, abs(DLavg), DLavg.real, \
      DLavg.imag ) )
  fin.close()
  fout.close()


