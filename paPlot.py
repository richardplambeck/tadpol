# paPlot.py
#
# compute polarization PA and percent vs time

import math
import time
import numpy
import subprocess
import shlex
import string


# Compute [PA,poli,sigmaPA,sigmaPoli] from [I,Q,U,sigmaI,sigmaQ,sigmaU]
def paCalc(I, Q, U, sigmaI, sigmaQ, sigmaU) :
  pa = (180./math.pi) * 0.5 * math.atan2(U,Q)
  poli = math.sqrt(U*U + Q*Q)
  dpolisq = (Q*Q*sigmaQ*sigmaQ + U*U*sigmaU*sigmaU)/(Q*Q + U*U)
  sigmaPoli = math.sqrt(dpolisq)
  dpasq = (0.25/pow(Q*Q+U*U,2.))
  dpasq = dpasq * (Q*Q*sigmaU*sigmaU + U*U*sigmaQ*sigmaQ) 
  sigmaPa = (180./math.pi) * math.sqrt(dpasq)
  return [pa, poli, sigmaPa, sigmaPoli]
  
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
     
# Run uvflux on selected data, return I, POLI, PA, and error estimates #
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

# --- this is the main routine that generates the table --- #
#  blocks of data are selected by source and time; extra is a string with any additional criteria
#   .e.g., extra="uvrange(20,100)"
# deltaMinutes is length of each time block, nreps=number of repeats for that time chunk
#    e.g., if source was observed for 12 minutes each time, use deltaMinutes=2, nreps=6

def makeTable( visFile='wide.av', \
               srcName='1733-130', \
               nint=10000, \
               extra="", \
               lineString='chan,1,1,8', \
               outFile='table.dat' ) :

  selectList = makeSelectList( visFile, srcName, nint, 2 )
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
 

# --- this is for vlbi processing - avg entire scan --- #
def makeTable2( visFile='wide.av', \
               srcName='1733-130', \
               extra="", \
               lineString='chan,1,1,8', \
               outFile='table.dat' ) :

  selectString="source(%s)" % srcName 
  if len(extra) > 0 :
    selectString = selectString + "," + extra 
  [I, sigmaI, pa, sigmaPa, poli, sigmaPoli, V, sigmaV] = getStokes( visFile, selectString, lineString )
  if (I > 0.) :		# I == 0 if uvflux returned an error
    polPercent = poli/I
    sigmaPolPercent = polPercent * math.sqrt( pow(sigmaPoli/poli,2) + pow(sigmaI/I,2) )
    fout = open(outFile, 'a')
    fout.write("  %8s  %8.3f %6.3f %8.3f %6.3f %7.1f %5.1f   %s\n" % \
      (srcName,I,sigmaI,polPercent,sigmaPolPercent,pa,sigmaPa,selectString) )
    fout.close()
 
# --- give Farhad limited data ---
def rework( infile, outfile ) :
   fin = open( infile, "r" )
   fout = open( outfile, "w" )
   for line in fin :
     if line.startswith("#") :
       fout.write( "%s" % line )
     else :
       fout.write( "%s %s %s" % ( line[0:10], line[18:43], line[89:] ) )
   fin.close()
   fout.close()

# retrieve the average UT, parallactic angle, and HA for time range specified by selectString
def getUTPAHA( visFile, selectString ) :
  p= subprocess.Popen( ( shlex.split('uvlist vis=%s select=%s,pol(LL) options=list recnum=1000000' % \
     ( visFile, selectString ) ) ), \
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
  if (navg > 0) :
    return [sumUT/float(navg), sumPA/float(navg), sumHA/float(navg)]
  else :
    return [ 0., 0., 0. ]

# make list of time ranges for a given source
# each time range contains up to nint integrations, but can contain up to 10% less if that
#   avoids a big time gap; this is designed to avoid weird groupings that happen because
#   of system temp measurements
# note: using select=time(t1,t2), a record is selected only if t >= t1 and t < t2;
#    thus, always add 1 sec to time of last record to make sure it will be selected

def makeSelectList( visFile, srcName, nint, maxGapMinutes ) :
  date1 = ""
  selectList = []
  scanList = []	    # string, e.g. '13MAR22:02:11:10.0'
  scanTime = []     # decimal minutes since 0 UT on date1

  # Create list of scan times for this source
  p= subprocess.Popen( ( shlex.split('uvindex vis=%s interval=0.01' %  visFile  ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  uvindexOutput = p.communicate()[0]
  lines = uvindexOutput.split("\n")
  for line in lines: 
    a = line.split()
    if (len(a) > 1) and (a[1] == srcName.upper() ) :
      scanList.append(a[0])
      [date,strhr,strmin,strsec] = a[0].split(":")
      hr = int(strhr)
      min = int(strmin)
      sec = float(strsec)
      decMinutes = 60.*hr + min + sec/60.
      if date1 == "" :
        date1 = date
      if (date != date1) :
        decMinutes = decMinutes + 1440.
      scanTime.append(decMinutes)

  # Parse the list
  nscans = 0
  for n in range(0,len(scanList)) :
    nscans = nscans + 1
    print "%s  %3d" % (scanList[n], nscans)
    if nscans == 1 :       
      t1 = scanList[n]               # set t1 if scan n is first scan in a group
    if (n < len(scanList)-1) :
      tgapnext = scanTime[n+1] - scanTime[n]    # gap before next scan
    # this is the last scan in a group if:
    #   ... quota is filled
    #   ... or this is the very last scan
    #   ... or quota is nearly filled and this scan is followed by big time gap
    if (nscans >= nint) or (n == len(scanList)-1) or (tgapnext > maxGapMinutes) :
      [date,strhr,strmin,strsec] = scanList[n].split(":")
      hr = int(strhr)
      min = int(strmin)
      sec = float(strsec) + 1.
      if sec >= 60. :
        sec = sec - 60.
        min = min + 1
      if min >= 60 :
        min = min - 60
        hr = hr + 1
      if (hr > 23) :
        print "### WARNING: WRAPPED THROUGH NEXT DAY ###"
        return [ ]
      t2 = "%s:%0.2d:%0.2d:%04.1f" % (date,hr,min,sec)
      oneString = "source(%s),time(%s,%s)" % (srcName,t1,t2)
      selectList.append( oneString ) 
      print oneString
      print " "
      nscans = 0
  return selectList

def RMcalc( fGHzLsb, PAlsb, fGHzUsb, PAusb ) :
  lambda_lsb_meters = 0.30/fGHzLsb
  lambda_usb_meters = 0.30/fGHzUsb
  deltaPhiRadians = (PAlsb - PAusb) * math.pi/180.
  RM = deltaPhiRadians/(pow(lambda_lsb_meters,2.) - pow(lambda_usb_meters,2.))
  print "RM = %.2e radians/m^2" % RM
  return RM

# read standard PA files for lsb, usb, compute pa(usb)-pa(lsb) and uncertainty
def computeRM( lsbfile, usbfile, rmfile ) :
  selectString = []
  lsbPA = []
  lsbSigma = []
  fin = open( lsbfile, "r" )
  for line in fin :
    a = line.split()
    if (not line.startswith("#")) and (len(a) == 12) :
      selectString.append( a[11] )
      lsbPA.append( float(a[7]) )
      lsbSigma.append( float(a[8]) )
  fin.close()    
  fin = open( usbfile, "r" )
  fout = open( rmfile, "w" )
  for line in fin :
    a = line.split()
    if (not line.startswith("#")) and (len(a) == 12) :
      for n in range (0, len(selectString) ) :
        if selectString[n] == a[11] :
          usbPA = float(a[7])
          usbSigma = float(a[8])
          dif = usbPA - lsbPA[n]
          sigma = math.sqrt(usbSigma*usbSigma + lsbSigma[n]*lsbSigma[n])
          RM = RMcalc( 217.75, lsbPA[n], 232.25, usbPA )
          print "%8s  %6s  %6s  %7.1f  %7.1f  %7.1f  %5.1f  %.2e" % \
            (a[0],a[1],a[2],lsbPA[n],usbPA,dif,sigma,RM)
          fout.write("%8s  %6s  %6s  %7.1f  %7.1f  %7.1f %5.1f  %.2e\n" % \
            (a[0],a[1],a[2],lsbPA[n],usbPA,dif,sigma,RM))
  fin.close()    
  fout.close()

def save() :
  makeTable( visFile='wide.av', srcName='SgrA', nint=3, extra="uvrange(20,1000)", lineString='chan,1,1,8', outFile='SgrA.nocal.30sec')
  makeTable( visFile='wide.av', srcName='1743-038', nint=3, extra="uvrange(20,1000)", lineString='chan,1,1,8', outFile='1743-038.nocal.30sec' ) 

def doit() :
  makeTable( visFile='wide.cal', srcName='SgrA', nint=3, extra="uvrange(20,1000)", lineString='chan,1,1,8', outFile='SgrA.30sec.2')
  makeTable( visFile='wide.cal', srcName='1743-038', nint=3, extra="uvrange(20,1000)", lineString='chan,1,1,8', outFile='1743-038.30sec.2' ) 

