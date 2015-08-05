# leakSolve.py 
# 27jun2012 - begin adapting from paPlot.py
# 05may2013 - deleted unused routines paCalc, parseLine, getStokes; as far as I can 
#   see, this code has nothing in common with paPlot.py
# 16jun2013 - delete code making separate leak files for each antenna - put all leakages into one Lk file
# 18oct2013 - add avgchan=0 option to average all channels in each window (this will be the default)
# 06nov2013 - change vis["avgchan"] to be a string, allowing "DSB", "2SB", or a number

import math
import time
import numpy
import subprocess
import shlex
import string
import sys
import struct

# this is a collection of routines (not objects!) that solve for or restore leakages;
#    leakages are written to or read from a disk file
# 
# begin with a fully calibrated data set; it may include many sources, and may have a
#    gains (and conceivably) passband items in it
#
# dictionary file controls operation - generally this will be created in calling routine
# vis1 = { "fileName"   : "wide.cal",                    (no default) (this should be "visName")
#          "selectStr"  : "source(3c279),-ant(7)",       (no default)
#          "LkName"     : "Lk.22mar2013",                (no default)
#          "avgchan"    : "5",                           (default = 0 -> avg across each corr section)
#          "interval"   : 5,                              
#          "refant"     : 8,                           
#          "optionStr"  : "circular,noxy,nopass,qusolve",
#          "flux"       : "10.",
#          "legend"     : "vlbi/22mar2013/3c279"          
#
# leakSolve.makeLk( vis )                                       - generates Lkfile (erase previous)
# leakSolve.restoreLk( LkFile, lineStr, outfile )               - restores leakage solution to outfile 
# leakSolve.stripout( Lkfile, antList, fstart, fstop, outfile ) - phased leakage sum for vlbi


def fillDefaults( vis ) :
  """Enter defaults for any missing dictionary items"""
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
  if not "avgchan" in vis :
    vis[ "avgchan" ] = "DSB"
  if not "legend" in vis :
    vis[ "legend" ] = ""
    
def delHd( fileName ) :
  """Delete existing leakage file"""
  p= subprocess.Popen( ( shlex.split('delhd in=%s/leakage' % ( fileName ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  
# copy (do not average) selected data from visFile into 'sstmp', with options=nopol
def copytoSStmp( visFile, selectStr ) :
  """copy data specified by selectStr from from visFile to sstmp; apply gain cal"""
  p= subprocess.Popen( ( shlex.split('rm -rf sstmp' ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  time.sleep(1)
  p=subprocess.Popen( ( shlex.split('uvcat vis=%s select=%s out=sstmp options=nopol' \
     % (visFile, selectStr) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  result = p.communicate()[0]
  print result

def makeFtable( vis ) :
  """Generates dictionary items "nstart", "nchan" (per window), "chfreq", "chwidth" (per chan)"""
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

def dumpGainRatio( visFile ) :
  '''Compute avg R/L gain ratio for each antenna in file vis'''
  print "dump gain ratio from file %s" % visFile
  g = []                 # accumulate gains in 1-D list of floats
  n = 0                  # number of time intervals
  RLratio = numpy.zeros(15)
  p= subprocess.Popen( ( shlex.split('gpplt vis=%s log=tmpGain' % visFile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  print result
  fin = open( 'tmpGain', 'r' )
  for line in fin :
    if (len(line) > 0) and (not line.startswith("#")) :
      a = line.split()
      if len(a) == 8 :  # begins new record
        n = n + 1       # increment line number
        del a[0:2]      # discard daynumber and UT
      for i in range(0,len(a)) :
        g.append( float( a[i] ) )
  gains = numpy.average( numpy.reshape( numpy.array(g), [n,-1] ), axis=0 )
  gRL = numpy.reshape( gains, [-1, 2] )
  for i in range(0,15) :
    if ( gRL[i][1] > 0. ) :
      RLratio[i] = gRL[i][0]/gRL[i][1]
  print RLratio 
  return RLratio
  
# --- particularly useful if averaging over multiple correlator windows --- #
def fminfmax( vis, ch1, ch2 ) :
  """find min and max freqs in channel range (ch1,ch2); uses vis[chfreq] and vis[chwidth] arrays"""
  fmin = 1000.   # GHz
  fmax = 0.      # GHz
  for nch in range(ch1,ch2+1) :
    f1 = vis["chfreq"][nch-1] - 0.5*vis["chwidth"][nch-1]
    f2 = vis["chfreq"][nch-1] + 0.5*vis["chwidth"][nch-1]
    if f1 < fmin : fmin = f1
    if f2 < fmin : fmin = f2
    if f1 > fmax  : fmax = f1
    if f2 > fmax  : fmax = f2
  return [fmin,fmax]

def addLk( vis, ch1, navg ) :
  """Solve for leakages for a line=chan,ch1,navg, append them to Lkfile"""
  lineString = "chan,1,%d,%d" % ( ch1, navg )
  ch2 = ch1 + navg - 1
  [fstart,fstop] = fminfmax( vis, ch1, ch2 )
  print "%s   %8.3f %8.3f" % (lineString, fstart, fstop)
  [ percentQ, percentU, DR, DL ] = runGpcal( vis, lineString=lineString )
  RLgainRatio = dumpGainRatio( vis["fileName"] )
  filename = vis["LkFile"]
  fout = open( filename, "a" )
  fout.write("#\n")
  for nant in range(0,15) :
    fout.write("C%02d %8.3f %8.3f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f   %s   %5.3f\n" % \
       ( nant+1, fstart, fstop, DR[nant].real, DR[nant].imag, DL[nant].real, \
       DL[nant].imag, percentQ, percentU, lineString, RLgainRatio[nant]) )
  fout.close()

def makeLk( visin ) :
  """Creates new Lkfile following prescription given in vis dictionary"""
  vis = visin.copy()
	# to make sure we don't mess up entries in vis, create temporary copy
  fillDefaults( vis )
  filename = vis["LkFile"]
  fout = open( filename, "w" )    # wipes out previous file!
  fout.write("# input data : %s\n" % vis["fileName"] )
  fout.write("# legend     : %s\n" % vis["legend"] )
  fout.write("# selectStr  : %s\n" % vis["selectStr"] )
  fout.write("# optionStr  : %s\n" % vis["optionStr"] )
  fout.write("# flux       : %s\n" % vis["flux"] )
  fout.write("# refant     : %d\n" % vis["refant"] )
  fout.write("# avgchan    : %s\n" % vis["avgchan"] )
  fout.write("# interval   : %.2f\n" % vis["interval"] )
  fout.close()
  copytoSStmp( vis["fileName"], vis["selectStr"] ) 
    # copies requested data to file sstmp, applying gain correction if present 
  vis["fileName"] = "sstmp"
  makeFtable( vis )	   # fill in nstart, nchan, chfreq, chwidth dictionary items
  nchlist = numpy.array( vis["nchan"] )   # needed only for DSB or 2SB
  nchtot = numpy.sum(nchlist)             # needed only for DSB or 2SB
  if vis["avgchan"].isdigit() :
    print "processing by correlator window" 
    for nstart,nchan in zip( vis["nstart"],vis["nchan"] ) :
      if nchan > 0 :      # skip corr sections with 0 chans
        avgchan = int(vis["avgchan"])
        if (avgchan > nchan) or (avgchan <= 0) :
          avgchan = nchan
        nextra = nchan - (nchan/avgchan)*avgchan  # number of leftover chans
        for n1 in range( nstart + nextra/2, nstart+nchan-avgchan+1, avgchan ) :
          print n1, n1+avgchan - 1
          delHd( vis["fileName"] ) 
          addLk( vis, n1, avgchan )
  elif vis["avgchan"] == "DSB" :
    print "computing DSB leakage"
    delHd( vis["fileName"] )
    addLk( vis, 1, nchtot )    # average all chans
  elif vis["avgchan"] == "LSB" :
    print "computing LSB eakage"
    delHd( vis["fileName"] )
    addLk( vis, 1, nchtot/2 )             # LSB
  elif vis["avgchan"] == "USB" :
    print "computing USB leakage"
    delHd( vis["fileName"] )
    addLk( vis, nchtot/2+1, nchtot/2 )    # USB
  elif vis["avgchan"] == "2SB" :
    print "computing LSB and USB leakages"
    delHd( vis["fileName"] )
    addLk( vis, 1, nchtot/2 )             # LSB
    delHd( vis["fileName"] )
    addLk( vis, nchtot/2+1, nchtot/2 )    # USB
  else :
    print "ERROR - vis[avgchan] = %s is unrecognized value" % vis["avgchan"]
   

def restoreLk( LkFile, lineStr, outfile ) :
  """Copies Lk specified by lineStr from LkFile to outfile"""
  print "restoreLk: copy leakages for lineStr = %s from %s to %s" % (lineStr, LkFile, outfile)
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
    leakStr = struct.pack( ">ii", 7, 0 )   # this is the preamble for all leakage items
    for i in range(0,60) :
      leakStr = leakStr + struct.pack( ">f", lk[i] )
    fout = open( outfile+"/leakage", "wb" )
    fout.write( leakStr )
    fout.close()

    # ---- writeleak no longer works, so had to replace this with code above ---- #
    # leakStr = ""
    # for n in range(0,60) :
    #   if n == 0 :
    #     leakStr = "%.3f" % lk[0]
    #   else :
    #     leakStr = leakStr + ",%.3f" % (lk[n])
    #print "leakStr=%s" % leakStr 
    # p= subprocess.Popen( ( shlex.split('/fringe2/plambeck/pol/tadpol/writeleak vis=%s leak=%s' % \
    # ( outfile, leakStr ) ) ), stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
    # result = p.communicate()[0]
    # if (len(result) > 1) : print result 
    # -------------------------------------------------------------------------- #

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

# --- plot correlator setups for Miriad data files listed in visList ---
def visualize( visList ) :
  fin = open( visList, "r" )
  for line in visList :
    if not line.beginswith("#") :
      try:
        vis["fileName"] = line
        makeFtable( vis )
        print vis["chfreq"] 
      except :
        print "failed to process %s" % vis["fileName"]
        
# --- plot R/L gain ratio --- #
def gainRatio( visFile, selectStr, refant=8, interval=10000 ) :
  parray = numpy.array( [2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.] )
  gR = numpy.zeros( 15 )
  gL = numpy.zeros( 15 )
  # use mfcal to solve for gains
  p= subprocess.Popen( ( shlex.split('mfcal vis=%s select=%s refant=%d interval=%0.1f' % \
     (visFile, selectStr, refant, interval ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  time.sleep(1)
  result = p.communicate()[0]
  print result
  # use gpplt to list the gains
  p= subprocess.Popen( ( shlex.split('gpplt vis=%s log=gainlog' % visFile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  time.sleep(1)
  result = p.communicate()[0]
  print result
  # read the gains file, and plot the gain ratios
  fin = open("gainlog", "r")
  ant = 1
  for line in fin :
    if not line.startswith("#") and len(line) > 1 :
      a = line.split()
      noffset = 0
      if len(a) == 8 :	   # begins new time
        ant = 1
        noffset = 2
      if (ant < 15) :
        for n in range(0,5,2) :
          gR[ant-1] = a[n+noffset]
          gL[ant-1] = a[n+1+noffset]
          ant = ant + 1
  print numpy.power(gL,parray)/numpy.power(gR,parray)
