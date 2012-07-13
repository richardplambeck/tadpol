# paPlot.py
#
# compute polarization PA and percent vs time

import math
import time
import numpy
import subprocess
import shlex
import string
import os

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

# --- specialized routine to measure beam polarizations in 24jun2012 Saturn dataset --- #
def offsetPAs() :
  offsets = [ [ '3C279',   '0.,0.'  ],      \
             [ '3C279O1', '-36.26,-.02' ], \
             [ '3C279O2', '-18.36,-.02' ], \
             [ '3C279O3', '-0.45,-9.02' ], \
             [ '3C279O4', '-0.45,8.98'  ], \
             [ '3C279O5', '17.46,-.02'  ], \
             [ '3C279O6', '35.36,-.02'  ] ]
  fout = open("beampol.dat", "a")
  fout.write("\n    name       offset      I    rmsI      PA  rmsPA    M/I      V/I\n")
  fout.close()
  for offset in offsets: 
    print offset[0]
    os.system("selfcal vis=lsb.cal select='source(%s)' offset=%s interval=0.1 refant=8" % (offset[0],offset[1]))
    os.system("selfcal vis=usb.cal select='source(%s)' offset=%s interval=0.1 refant=8" % (offset[0],offset[1]))
    p= subprocess.Popen( ( shlex.split("uvflux vis=lsb.cal,usb.cal select='source(%s)' line=chan,1,1,96 offset=%s stokes=I,Q,U,V" \
       % (offset[0],offset[1]) ) ), \
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
    fout = open("beampol.dat", "a")
    fout.write("%8s %12s    %5.2f  %4.2f  %7.2f  %4.2f   %5.3f   %6.3f\n" % ( offset[0],offset[1], I, rmsI, pa, sigmaPa, poli/I, V/I ) )
    fout.close()

# --- this is the main routine that generates the table --- #
#  blocks of data are selected by source and time; extra is a string with any additional criteria
#   .e.g., extra="uvrange(20,100)"
# deltaMinutes is length of each time block, nreps=number of repeats for that time chunk
#    e.g., if source was observed for 12 minutes each time, use deltaMinutes=2, nreps=6

def makeTable( visFile='wide.av', srcName='1733-130', extra="", deltaMinutes=3., nreps=1, lineString='chan,1,1,8', outFile='table.dat' ) :
  p= subprocess.Popen( ( shlex.split('uvindex vis=%s' % visFile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  uvindexOutput = p.communicate()[0]
  [selectList,timeList] = parseIndex( uvindexOutput, srcName, deltaMinutes, nreps )
  print "\n"
  print selectList
  fout = open(outFile, 'a')
  fout.write("#\n# ------ visFile = %s, lineString = %s ------ #\n" % (visFile,lineString))
  fout.write("#  dechr   parang     HA        S    sigma     poli  sigma     PA  sigma      V    sigma   selectString\n")
  fout.close()
  n = -1 
  for selectString in selectList :
    if len(extra) > 0 :
      selectString = selectString + "," + extra 
    print " "
    print selectString 
    n = n + 1
    [I, sigmaI, pa, sigmaPa, poli, sigmaPoli, V, sigmaV] = getStokes( visFile, selectString, lineString )
    [UT,PA,HA] = getUTPAHA( visFile, selectString ) 
    if (I > 0.) :		# I == 0 if uvflux returned an error
      fout = open(outFile, 'a')
      fout.write("%8.3f %8.2f %8.3f %8.3f %6.3f %8.3f %6.3f %7.1f %5.1f %8.3f %6.3f   %s\n" % \
	    (UT,PA,HA,I,sigmaI,poli,sigmaPoli,pa,sigmaPa,V,sigmaV,selectString) )
      fout.close()
 
# retrieve the average UT, parallactic angle, and HA for time range specified by selectString
#
def getUTPAHA( visFile, selectString ) :
  p= subprocess.Popen( ( shlex.split('uvlist vis=%s select=%s options=list recnum=10000' % ( visFile, selectString ) ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  uvlistOutput = p.communicate()[0]
  navg = 0
  sumUT = 0.
  sumHA = 0.
  sumPA = 0.
  lines = uvlistOutput.split("\n")
  for line in lines:
    a = line.split()
    if (len(a) == 19) :
      navg = navg + 1
      sumUT = sumUT + float(a[2])
      sumHA = sumHA + float(a[4])
      sumPA = sumPA + float(a[16])
  return [sumUT/float(navg), sumPA/float(navg), sumHA/float(navg)]

def test1() :
  getPAHA( 'bllac.lsb', 'source(BLLAC),time(12MAY07:09:42:14.5,12MAY07:09:43:15.5),ant(1)(2)' )

# generate table of selectStrings from indexFile; each contains source and time range
#   e.g., 'source(BLLAC),time(12MAY07:09:42:14.5,12MAY07:09:43:15.5)'
# unfortunately uvindex.log contains only the starting time of each observation, not 
#   the total duration; so, this routine generates nreps selectStrings
def parseIndex( uvindexOutput, srcName, deltaMinutes, nreps ) :
  selectList = []
  timeList = []
  #lines = uvindexOutput.split("\n")
  try :
    fin = open( "uvindex.log", "r")
  except :
    print "ERROR: uvindex.log does not exist"
    exit() 
  for line in fin:
    print line
    if ( (not line.startswith("#")) and (string.count(line, ":") == 3) ) :
      print line
      a = line.split()
      if a[1] == srcName :
        t1 = a[0]
        for n in range(0,nreps) :
          [t2,t3] = trange(t1,n,deltaMinutes)
          oneString = "source(%s),time(%s,%s)" % (srcName,t2,t3)
          selectList.append( oneString )
          timeList.append( dechrs(t2) + 0.5*deltaMinutes/60. )
    if (string.count(line,"-") > 20) :		# signifies end of relevant section of uvindex output
	  break
  return selectList,timeList

# convert timeString to decimal hrs, without taking into account the day
def dechrs( timeString ) :
  [d1,h1,m1,s1] = timeString.split(":")
  return int(h1) + float(m1)/60. + float(s1)/3600.
          
# generates start, stop time strings (e.g., "12MAY04:12:01:44.0")
# these specify a time interval of deltaMinutes, beginning at
#   startTimeString + nincr*deltaMinutes
def trange( startTimeString, nincr, deltaMinutes ) :
  [d1,h1,m1,s1] = startTimeString.split(":")
  s2 = float(s1) + nincr*60.*deltaMinutes
  m2 = int(m1)
  h2 = int(h1)
  while (s2 >= 60.) :
    s2 = s2 - 60.
    m2 = m2 + 1
  while (m2 >= 60) :
    m2 = m2 - 60
    h2 = h2 + 1
  if (h2 > 23) :
    print " error - goes through 1 day "

  s3 = float(s1) + (nincr+1)*60.*deltaMinutes + 1.
  m3 = int(m1)
  h3 = int(h1)
  while (s3 >= 60.) :
    s3 = s3 - 60.
    m3 = m3 + 1
  while (m3 >= 60) :
    m3 = m3 - 60
    h3 = h3 + 1
  if (h3 > 23) :
    print " error - goes through 1 day "
 
  t2 = "%s:%0.2d:%0.2d:%03.1f" % (d1,h2,m2,s2)
  t3 = "%s:%0.2d:%0.2d:%03.1f" % (d1,h3,m3,s3)
  return [t2,t3]

def savethese() :
  makeTable( visFile='wide.lsb', srcName='SGRA', extra='-ant(12),uvrange(20,1000)', deltaMinutes=1., nreps=12, lineString='chan,1,1,4', outFile='sgra.1min.lsb' ) 
  makeTable( visFile='wide.usb', srcName='SGRA', extra='-ant(12),uvrange(20,1000)', deltaMinutes=1., nreps=12, lineString='chan,1,1,4', outFile='sgra.1min.usb' ) 
  makeTable( visFile='wide.lsb', srcName='3C279', extra='-ant(12)', deltaMinutes=1., nreps=5, lineString='chan,1,1,4', outFile='3c279.lsb' ) 
  makeTable( visFile='wide.usb', srcName='3C279', extra='-ant(12)', deltaMinutes=1., nreps=5, lineString='chan,1,1,4', outFile='3c279.usb' ) 
  makeTable( visFile='wide.lsb', srcName='SGRA', extra='-ant(12),uvrange(20,1000)', deltaMinutes=4, nreps=3, lineString='chan,1,1,4', outFile='sgra.4min.lsb' ) 
  makeTable( visFile='wide.usb', srcName='SGRA', extra='-ant(12),uvrange(20,1000)', deltaMinutes=4, nreps=3, lineString='chan,1,1,4', outFile='sgra.4min.usb' ) 

def doit() :
  makeTable( visFile='wide.lsb', srcName='1924-292', extra='-ant(12)', deltaMinutes=3., nreps=1, lineString='chan,1,1,4', outFile='1924.lsb' ) 
  makeTable( visFile='wide.usb', srcName='1924-292', extra='-ant(12)', deltaMinutes=3., nreps=1, lineString='chan,1,1,4', outFile='1924.usb' ) 
  makeTable( visFile='wide.lsb', srcName='1733-130', extra='-ant(12)', deltaMinutes=3., nreps=1, lineString='chan,1,1,4', outFile='1733.lsb' ) 
  makeTable( visFile='wide.usb', srcName='1733-130', extra='-ant(12)', deltaMinutes=3., nreps=1, lineString='chan,1,1,4', outFile='1733.usb' ) 
  makeTable( visFile='wide.lsb', srcName='3C286', extra='-ant(12)', deltaMinutes=3., nreps=5, lineString='chan,1,1,4', outFile='3C286.lsb' ) 
  makeTable( visFile='wide.usb', srcName='3C286', extra='-ant(12)', deltaMinutes=3., nreps=5, lineString='chan,1,1,4', outFile='3C286.usb' ) 

# --- generates wip commands to draw horiz bars for each corr window on a frequency plot --- #
def spectraPlot( infile, outfile ) :
  nchan = []
  fstart = []
  finterval = []
  p= subprocess.Popen( ( shlex.split('uvlist vis=%s options=spectra' % infile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  lines = result.split("\n")
  for line in lines :
    a = line.split()
    if line.startswith("number of channels") :
      for n in range (4, len(a)) :
        nchan.append( int(a[n]) )
    if line.startswith("starting frequency") :
      for n in range (3, len(a)) :
        fstart.append( float(a[n]) )
    if line.startswith("frequency interval") :
      for n in range (3, len(a)) :
        finterval.append( float(a[n]) )
  fout = open( outfile, "w" )
  fout.write("mtext t -2 .05 0 %s\n" % infile )
  for n in range(0, len(nchan)) :
    if nchan[n] > 0 : 
      fstop = (nchan[n]-1)*finterval[n] + fstart[n]
      fwidth = abs( fstop - fstart[n] ) 
      print "%3d  %7.3f  %7.3f" % (nchan[n], fstart[n], fstop)
      if (fwidth < .4) :
        fout.write("color 2\n")
      fout.write("move %7.3f 0.\n" % fstart[n] )
      fout.write("lwidth 12\n")
      fout.write("draw %7.3f 0.\n" %  fstop )
      fout.write("lwidth 1\n")
      fout.write("move %7.3f 0.002\n" % ( (fstart[n] + fstop)/2. ) )
      fout.write("putlabel 0.5 %d\n" % (n+1) )
      if (fwidth < .4) :
        fout.write("color 1\n")
  fout.close() 
    

  

