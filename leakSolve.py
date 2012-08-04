# leakSolve.py 
# 27jun2012 - begin adapting from paPlot.py
#

import math
import time
import numpy
import subprocess
import shlex
import string


# --- input I, Q, U, sigmas; return polm, pa, sigmas --- #
# note that sigmaI, sigmaQ, sigmaU are uncertainties of the mean
def paCalc(I, Q, U, sigmaI, sigmaQ, sigmaU) :
  pa = (180./math.pi) * 0.5 * math.atan2(U,Q)
  poli = math.sqrt(U*U + Q*Q)
  dpolisq = (Q*Q*sigmaQ*sigmaQ + U*U*sigmaU*sigmaU)/(Q*Q + U*U)
  sigmaPoli = math.sqrt(dpolisq)
  dpasq = (0.25/pow(Q*Q+U*U,2.))
  dpasq = dpasq * (Q*Q*sigmaU*sigmaU + U*U*sigmaQ*sigmaQ) 
  sigmaPa = (180./math.pi) * math.sqrt(dpasq)
  return [pa, poli, sigmaPa, sigmaPoli]
  
# --- extracts Stokes amplitude and rms from one line of uvflux output --- #
def parseLine( line ) :
  print " "
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
     
# --- runs uvflux on selected data, returns I, POLI, PA, and error estimates --- #
def getStokes( infile, selectString, lineString ) :
  p= subprocess.Popen( ( shlex.split('uvflux vis=%s select=%s line=%s stokes=I,Q,U,V' \
     % (infile, selectString, lineString) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  lines = result.split("\n")
  I = rmsI = pa = sigmaPa = poli = sigmaPoli = V = rmsV = 0.	# in case we get error
  nlines = len(lines)
  if nlines >= 10 :
    [src, I, rmsI, ncorrs] = parseLine( lines[nlines-6] )
    [dummy, Q, rmsQ, ncorrs] = parseLine( lines[nlines-5] )
    [dummy, U, rmsU, ncorrs] = parseLine( lines[nlines-4] )
    [dummy, V, rmsV, ncorrs] = parseLine( lines[nlines-3] )
    [pa, poli, sigmaPa, sigmaPoli] = paCalc(I, Q, U, rmsI, rmsQ, rmsU) 
  return [ I, rmsI, pa, sigmaPa, poli, sigmaPoli, V, rmsV ]


def runGpcal( visFile='', \
              selectString='', \
              lineString='', \
              fluxString='', \
              refant=1, \
              interval=0.5, \
              optionString='circular,qusolve,noxy,nopass') :
  DR = numpy.zeros( 15, dtype=complex)
  DL = numpy.zeros( 15, dtype=complex)
  if fluxString == "" :
    fluxString = "1."
  p= subprocess.Popen( ( shlex.split('gpcal vis=%s select=%s line=%s flux=%s refant=%d interval=%0.1f options=%s' % \
     (visFile, selectString, lineString, fluxString, refant, interval, optionString ) ) ), \
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
      print line
      ant = int(line[4:6])
      if ant < 16 :
        DR[ant-1] = float(line[16:22]) + 1.j * float(line[23:29])
        DL[ant-1] = float(line[32:38]) + 1.j * float(line[39:45])
  return [percentQ,percentU,DR,DL]
       
def delHd( visFile ) :
  p= subprocess.Popen( ( shlex.split('delhd in=%s/leakage' % (visFile) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  
def makeFtable( visFile ) :
  nstart = []
  nchan = []
  fstart = []
  fstep = []
  chfreq = []
  chwidth = []
  p= subprocess.Popen( ( shlex.split('uvlist options=spectra vis=%s' % visFile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  uvlistOutput = p.communicate()[0]
  lines = uvlistOutput.split("\n")
  for line in lines :
    #print line
    a = line.split()
    if "starting channel" in line :
      for n in range( 3, len(a) ) :
        nstart.append( int(a[n]) ) 
    if "number of channels" in line :
      for n in range( 4, len(a) ) :
        nchan.append( int(a[n]) )
    if "starting frequency" in line :
      for n in range( 3, len(a) ) :
        fstart.append( float(a[n]) ) 
    if "frequency interval" in line :
      for n in range( 3, len(a) ) :
        fstep.append( float(a[n]) )
  #print "# ========================================================================= #"
  #print nstart, nchan, fstart, fstep      
  #print "# ========================================================================= #"
  for n in range(0, len(nstart) ) :
    if nchan[n] > 0 :
      for i in range( 0, nchan[n] ) :
         chfreq.append(fstart[n] + i*fstep[n])
         chwidth.append(fstep[n])
  #for ichn in range(0, len(chfreq) ) :
  #  print ichn+1, chfreq[ichn], chwidth[ichn]
  return [chfreq,chwidth]

def addLeak( visFile, ch1, nchan ) :
  [chfreq,chwidth] = makeFtable( visFile )
  lineString = "chan,1,%d,%d" % (ch1,nchan)
  ch2 = ch1 + nchan - 1
  fstart = chfreq[ch1-1] - 0.5*chwidth[ch1-1]
  fstop = chfreq[ch2-1] + 0.5*chwidth[ch2-1]
  print "%s   %8.3f %8.3f" % (lineString, fstart, fstop)
  [percentQ,percentU,DR,DL] = runGpcal( visFile=visFile, \
              selectString='', \
              lineString=lineString, \
              fluxString='', \
              refant=1, \
              interval=5, \
              optionString='circular,qusolve,noxy,nopass') 
  for nant in range(0,15) :
    filename = "leak" + "%d" % (nant+1)
    fout = open( filename, "a" )
    fout.write("%8.3f %8.3f %8.3f %6.1f %8.3f %6.1f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f    %s\n" % \
       ( fstart, fstop, numpy.abs(DR[nant]), numpy.angle(DR[nant],deg=True), \
       numpy.abs(DL[nant]), numpy.angle(DL[nant],deg=True), DR[nant].real, \
       DR[nant].imag, DL[nant].real, DL[nant].imag, percentQ, percentU, lineString) )
    fout.close()


def makeTable( visFile, nchan ) :
  for nwin in range (0,8) :
    n1 = 24*nwin + 1
    while (n1+nchan-1) < 24*(nwin+1)+1 :
      delHd( visFile )
      addLeak( visFile, n1, nchan )
      n1 = n1 + nchan
  
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

    # plot smooth curve
      if not hist :
        if abs(prevx - fmean) < (1.1*fwidth) :
          fout.write("draw %8.3f %s\n" % (fmean,a[ncol]) )
        else :
          fout.write("move %8.3f %s\n" % (fmean,a[ncol]) )
        if dots :
          fout.write("dot\n")  
        prevx = fmean

    # plot histogram
      else :
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
  
def doit() :
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

    

  
  
