# vCal.py
#
# computes SEFD for phased and unphased antennas during vlbi run
# 2013 version

import math
import cmath
import time
import numpy
import subprocess
import shlex
import timePlots
import datetime

# --- ants phased in 2013
# antList = [2,4,5,6,8,9,13,14]
 
# --- ants phased in 2015
# antList=[2,3,4,5,6,13,14,15]

# --- return time and gainComplex arrays as read from gplist output --- #
def getGains( infile ) :
  time = []
  gain = [] 
  p= subprocess.Popen( ( shlex.split('gplist vis=%s options=all' % infile) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT) 
  result = p.communicate()[0]
  lines = result.split("\n")
  ngains = (len(lines) - 3)/23
    # caution: this presumes 23 antennas will be listed for each time, and that
    # gains output has 3 header lines (and one blank line at the end ?)
  gainComplex = numpy.zeros( (ngains,15), dtype=complex )
  ng = -1 
  for n in range(3, len(lines) ) :
    a = lines[n].split()
    if ( (len(a) > 0) and (a[0] != "Ant") ) : 
      ng = ng + 1
      time.append( a[0] )
      nant = int(a[2])
      if ( nant != 1 ) :
        print "getGains error - unexpected ant number"
      try :
        gainComplex[ng][nant-1] = float(a[5]) + 1j * float(a[6])
      except :
        gainComplex[ng][nant-1] = 0 + 1j * 0 
    elif ( len(a) > 0 ) :
      nant = int(a[1]) 
      if ( nant < 16 ) :
        gainComplex[ng][nant-1] = float(a[4]) + 1j * float(a[5])
  return [time, gainComplex]

# new on 5/15/2015
# smooth the gain amplitudes with running boxcar average
# I decided NOT to smooth the gains in this way because of rapid gain variations
#   on C14 caused by temperature fluctuations

def smoothGainAmps( time, gainComplex, hrs=1. ) :
  t = numpy.empty( len(time), dtype=float )
  for n in range (0, len(time) ) :
    t[n] = dechrs( time[n] )
    # ... create decimal time array 
  gAmps = numpy.abs( gainComplex )
    # ... create array of gain amplitudes
  source = getSource ( "source.log", t )
  for n in range (0, len(time)-1 ) :
    if source[n] == "SGRA" :
      gAmps[n] = 0.
      # ... set SGRA gains to 0
  gAmps1 = numpy.ma.masked_outside( gAmps, 0.1, 5 )
    # ... toss out any gain amplitudes < 0.1 or > 5.
    # ... since SGRA gains were 0, they will be masked off too
  smoothGain = numpy.empty( [len(time),15] )
  istart = istop = 0
  for n in range (0, len(time) ) :
    #istart = n - 20
    #if istart < 0 :
    #  istart = 0
    #istop = n + 20
    #if istop > (len(time)-1) :
    #  istop = len(time)-1
    tcenter = t[n]
    while t[istart] < (tcenter - hrs/2.) :
      istart = istart + 1
      if istart > len(t) :
        print "ERROR: istart off end of time array"
    while (t[istop] < tcenter + hrs/2.) and (istop < (len(t)-1) ) :
      istop = istop + 1 
    smoothGain[n] =  numpy.ma.median( gAmps1[istart:istop], axis=0 )  
    print time[n], numpy.array_str( smoothGain[n], precision=3, max_line_width=200 )
  # replace any crazy values with average for that antenna
  gAmps2 = numpy.ma.masked_outside( smoothGain, 0.2, 3. )
  gAmpsAvg = numpy.ma.median( gAmps2, axis=0 )
  for i in range(0,15) :
    if numpy.isnan( gAmpsAvg[i] ) : 
      gAmpsAvg[i] = 0.
      # ... to prevent error message for C10
  print "\ngAmpsAvg = ", numpy.array_str( gAmpsAvg, precision=3, max_line_width=200 )
  for n in range( 0, len(time) ) :
    for i in range(0,15) :
      if smoothGain[n][i] == 0. :
        smoothGain[n][i] = gAmpsAvg[i]
  print numpy.array_str( numpy.ma.median( gAmps1[istart:istop], axis=0), precision=3, max_line_width=200 )
  timePlots.plotGains( t, gAmps, smoothGain )
  
# new on 5/18/2015
# (1) rescale antenna gain amplitudes on SgrA to median gain amplitudes of neighboring data
# (2) set C1 gains during scans equal to median over the full night (or to override value)
# 
# note that gainComplex = 0+0j for data that were flagged for any reason (C1 during scans,
#   pointing bad, shadowed); these are left alone except for C1 during scans, for data
#   records where C1's pointing was OK 
# 
# ... time = ascii time strings
# ... gainComplex is [ntimes x 15] matrix of complex gains from selfcal
# ... (t1str, t2str) is time interval to use for computing gain amp for Sgr
# ... C1amp = override value of C1 gain amp to use if desired
#
def SgrGain( time, gainComplex, t1str, t2str, C1amp=None, plotName="gain.pdf", foutCalc=None ) :

    outGain = numpy.array(gainComplex)
      # ... default is to leave gain alone

  # create decimal time array
    t = numpy.empty( len(time), dtype=float )
    for n in range (0, len(time) ) :
      t[n] = dechrs( time[n] )
    source = getSource( "source.log", t )
      # returns len(t) list of sources
    axisrms = getVar30( "axisrms.log", t )
      # returns len(t) x 15 matrix of rms pointing errors

  # derive median gain over the whole night, omitting crazy values
    print "... derive median gain amps over the whole night, omitting outliers and SgrA"
    print "... derive SgrA gain amplitudes from medians over time range (%s,%s)" % (t1str,t2str)
    foutCalc.write("# derive SgrA gain amplitudes from medians over time range (%s,%s)\n" % (t1str,t2str))
    gAmps = numpy.abs( gainComplex )
      # ... create array of gain amplitudes
    for n in range( 0, len(time) ) :
      if source[n] == "SGRA" :
        gAmps[n] = numpy.zeros( 15 )
          # ... set SGRA gains to 0 so they will be masked and not used to compute median
    gAmpsMasked = numpy.ma.masked_outside( gAmps, 0.1, 5 )
      # ... toss out any gain amplitudes < 0.1 or > 5. (SGRA will be masked off too)
    medianAmp = numpy.ma.median( gAmpsMasked, axis=0 )
      # ... these are the MEDIAN gain amplitudes for the entire night, ignoring masked values
    if C1amp :
      print "... using override gain amp of  %.3f for C1, instead of computed median %.3f" % \
      foutCalc.write( "#\n# using override gain amp of  %.3f for C1, instead of computed median %.3f\n" % \
        C1amp, medianAmp[0] )
      medianAmp[0] = C1amp
    print "\nmedian gains:  ", numpy.array_str( medianAmp, precision=3, max_line_width=200 )
    foutCalc.write("# median gains:   %s\n" % numpy.array_str( medianAmp, precision=3, max_line_width=200 ) )

  # derive SgrA gains over the interval (t1str,t2str)
    t1 = dechrs( t1str )
    t2 = dechrs( t2str )
  # find array indices (istart,istop) that correspond to t1str and t2str
    istart = istop = 0
    for n in range (0, len(time) ) :
      if t[n] < t1 :
        istart = istart + 1
        istop = istop + 1 
      elif t[n] < t2 :
        istop = istop + 1
    SgrAmp =  numpy.ma.median( gAmpsMasked[istart:istop], axis=0 )  
    SgrAmp[0] = medianAmp[0]
      # for C1 we have very few gain values, so use medianAmp for Sgr as well
    print "SgrA gains:    ", numpy.array_str( SgrAmp, precision=3, max_line_width=200 )
    print " "
    foutCalc.write("# Sgr gains:      %s\n#\n" % numpy.array_str( medianAmp, precision=3, max_line_width=200 ) )

  # now rescale amplitudes of SgrA gains, fill in C1 gain during scans
    msk = numpy.ma.getmask( gAmpsMasked )
    for n in range( 0, len(time) ) :
      for i in range( 0, 15 ) :
        if msk[n][i] :
          gA = medianAmp[i]
          if (t[n] > t1) and (t[n] < t2) :
            gA = SgrAmp[i]
          if (i==0) and (axisrms[n][i] < 10.) :
            outGain[n][i] = gA + 0j
          else :
            gAraw = abs(gainComplex[n][i])
            if (numpy.isnan(gAraw)) or (gAraw == 0.) :
              outGain[n][i] = 0. + 0j
            else :
              outGain[n][i] = gA/gAraw * gainComplex[n][i]
          #if source[n] == "SGRA" :
          #  print n, i, source[n], gainComplex[n][i], outGain[n][i]

  # create pdf plot if desired
    if plotName :
      print "... generating plot %s" % plotName
      newAmps = numpy.abs( outGain )
      gAmps = numpy.abs( gainComplex )
       # ... recompute these since Sgr gAmps were set to zero for masking
      timePlots.plotGains( t, gAmps, newAmps, plotName )

    return outGain   

# --- read varplt log containing axisrms (duplicate values for each ant)
# ... t is the decimal time array
def getVar30( infile, t ) :
  varlist = []
  fin = open( infile, "r" )
  n1 = 0
  nline = 0
  nskip = 0
  for line in fin:
    if ( not line.startswith("#") ) :
      a = line.split()
      if ( len(a) == 6) and (nline < len(t)) : 
        # print a[1], t[nline]
        if ( (nline >= len(t)) or (abs(dechrs(a[1]) - t[nline]) > .0008)) :
          #print "skipping %s record at UT %s - time not in gains array" % (infile,a[1])
          nskip = nskip + 1
        else :
          nline = nline + 1
          n1 = 2
          varlist.append( a[2] )
          varlist.append( a[4] )
      elif ( (len(a) == 4) and (n1 < 13) ) :
        n1 = n1 + 2
        varlist.append( a[0] )
        varlist.append( a[2] )
      elif ( (len(a) == 4) and (n1 == 14) ) :
        n1 = n1 + 1
        varlist.append( a[0] )
  fin.close()
  vl = numpy.array( varlist, dtype=float )
  # print "... skipped %d input lines in %s because times were not in gains array" % (nskip, infile)
  return numpy.reshape( vl, (-1,15) )

# --- read varplt log containing tsys or elev, 6 lines per vis record, 4 values per line --- #
# ... t is the decimal time array
def getVar15( infile, t ) :
  varlist = []
  fin = open( infile, "r" )
  n1 = 0
  nline = 0
  nskip = 0
  for line in fin:
    if ( not line.startswith("#") ) :
      a = line.split()
      if ( len(a) == 6) : 
        # print a[1], t[nline]
        if nline >= len(t) :
          # print "skipping %s record at UT %s - time not in gains array" % (infile,a[1])
          nskip = nskip + 1 
        elif abs(dechrs(a[1]) - t[nline]) > .0008 :
          # print "skipping %s record at UT %s - time not in gains array" % (infile,a[1])
          nskip = nskip + 1 
        else :
          nline = nline + 1
          n1 = 4
          for n in range(2,6) :
            varlist.append( a[n] )
      elif ( (len(a) == 4) and (n1 < 12) ) :
        n1 = n1 + 4
        for n in range(0,4) :
          varlist.append( a[n] )
      elif ( (len(a) == 4) and (n1 == 12) ) :
        n1 = n1 + 4
        for n in range(0,3) :
          varlist.append( a[n] )
  fin.close()
  vl = numpy.array( varlist, dtype=float )
  # print "... skipped %d input lines in %s because times were not in gains array" % (nskip, infile)
  return numpy.reshape( vl, (-1,15) )

# --- read varplt log containing single variable like rmspath or tau225 --- "
def getVar( infile, t ) :
  var = []
  fin = open( infile, "r" )
  nline = 0
  nskip = 0
  for line in fin:
    if ( not line.startswith("#") ) :
      a = line.split()
      if ( nline >= len(t) ) or ( abs(dechrs(a[1]) - t[nline]) > .0008 ) :
        # print "skipping %s record at UT %s = %.5f; next gains time = %.5f" % (infile,a[1],dechrs(a[1]),t[nline])
        nskip = nskip + 1
      else :
        var.append( a[2] )
        nline = nline + 1
  fin.close()
  # print "... skipped %d input lines in %s because times were not in gains array" % (nskip, infile)
  return numpy.array( var, dtype=float )

# --- read source names from 'uvindex interval=0.1' log file --- #
def getSource( infile, t ) :
  nline = 0
  nl = 0
  source = []
  fin = open( infile, "r" )
  for line in fin:
    nl += 1
    if nl > 7  :
      a = line.split()
      if ( a[1] == "Total" ) :
        fin.close()
        return source
      else :
        date,hr,min,sec = a[0].split(":")
        dtime = float(hr) + float(min)/60. + float(sec)/3600.
        if nline < len(t) :
          if ( abs(dtime - t[nline]) > .0003) or (nline >= len(t)) :
            # print "skipping %s record at UT %s - time not in gains array" % (infile,a[0])
            continue 
          else :
            source.append( a[1] )
            nline = nline + 1

# --- read source names from 'uvindex interval=0.1' log file --- #
def getFlux( infile, t ) :
  nline = 0
  nl = 0
  flux = []
  fin = open( infile, "r" )
  for line in fin:
    nl += 1
    if nl > 10  :
      a = line.split()
      decHrs = 24. * float(a[0])
      if ( abs(decHrs - t[nline]) > .0003) or (nline >= len(t)) :
        # print "skipping %s record at UT %s - time not in gains array" % (infile,a[0])
        continue 
      else :
        flux.append( float(a[1]) )
        nline = nline + 1
  return flux

def efficiency( s, s_ideal ) :
  if (s > 0.) :
    return s_ideal/s
  else :
    return 0. 

# --- compute avg of variable over a scan --- #
def avgVar( t1, t2, t, var ) : 
  navg = 0
  avg = 0.
  for n in range( 0, len(t) ) :
    if ( (t[n] >= t1) and (t[n] <= t2) ) :
      avg += var[n]
      navg += 1
  if (navg > 0) :
    return avg/navg
  else :
    return 0.

# --- wrapper for avgSEFD that uses default gains if there are no pre- or post-scan data --- #
# ... this should no longer be necessary, since SgrGain interpolates all gains, but
# ...    I will keep using it just in case

def scanSEFD( foutCalc, foutRbyR, t1, t2, src, t, source, gain, defgain, tsys, antList ) :
  s = avgSEFD( foutCalc, t1, t2, src, t, source, gain, tsys, antList ) 
  if (s > 0.) :
    if foutRbyR :
      foutRbyR.write("%8.0f " % (100.*round(s/100.)) )
  else :
    foutCalc.write("\n\n**** RECALCULATING using default gains because no valid gains ****\n")
    s = avgSEFD( foutCalc, t1, t2, src, t, source, defgain, tsys, antList )  
    if foutRbyR :
      foutRbyR.write("%8.0f*" % (100.*round(s/100.)) )
  return s


# SEFD returns SEFD for a single integration, for 1 antenna or for a phased array
# ... onegain = 1 x 15 complex array
# ... onetsys = 1 x 15 real array
# ... antList = list of antennas
# no printed or written output
# antennas with gain of 0 are assumed to have been flagged bad because they were not
#   pointed at source, or not locked; they add noise, but not signal, to phased sum
# if all antennas in list have gain of 0, return very large SEFD = 9999999.
#
def SEFD( onegain, onetsys, antList ) :
  jyperkNominal = [65.,65.,65.,65.,65.,65.,145.34,145.34,145.34,145.34,145.34,145.34,145.34,145.34,145.34]
  sum = 0.+ 0j
  for ant in antList :
    magg = numpy.abs( onegain[ant-1] )
    if ( (magg > 0.) and (onetsys[ant-1] > 10. ) ) :
      sum += 1./( onegain[ant-1] * math.sqrt( jyperkNominal[ant-1] * onetsys[ant-1] ) )
	      # note/ that denominator is vector amplitude, not power
  if abs(sum) > 0. :
    return len(antList)/(numpy.abs(sum) * numpy.abs(sum))  
  else :
    return 9999999.

def testSEFD( glist, tsyslist, antlist ) :
  gnComplex = numpy.zeros( 15, dtype=complex )
  tsys = 300.*numpy.ones( 15, dtype=float )
  for gn, ts, ant in zip (glist,tsyslist,antlist) :
    gnComplex[ant-1] = gn[0] * cmath.exp(1j*math.radians(gn[1]))
    tsys[ant-1] = ts
  print " "
  print "glist = ", glist
  print "tsys = ", tsyslist
  print "antlist = ", antlist
  print "SEFD = ", SEFD( gnComplex, tsys, antlist )
    
# AscanSEFD returns list of SEFD values, one for each integration in time range t1-t2
# the final value in the list is the average for the entire scan
# also writes out detailed info to file foutCalc, in case someone wants to check the calculation

def AscanSEFD( foutCalc, t1, t2, title, t, gain, tsys, antList ) :
    foutCalc.write("#\n# %s  antList = %s\n" % (title, antList))
  # write header for foutCalc section
    foutCalc.write("# UThrs ")
    for ant in antList :
      if ant < 10 :
        foutCalc.write( "     --- C%d ---" % ant ) 
      else :
        foutCalc.write( "    --- C%2d ---" % ant ) 
    foutCalc.write("      SEFD\n")
    navg = 0
    avg = 0.
    sList = []

  # process each integration within the time interval 
    for n in range( 0, len(t) ) :
      #if ( (t[n] >= t1) and (t[n] <= t2) and (src == source[n]) ) :
      if (t[n] >= t1) and (t[n] <= t2) :
        s = SEFD( gain[n], tsys[n], antList ) 
          # ... note that gain[n] and tsys[n] are (1 x 15) arrays
        sList.append( s )

      # write 1 line with gains and scan SEFD --- #
        str = "%8.4f" % t[n]
        for ant in antList :
          str = str + "  %8.3f %4.0f" % ( abs(gain[n][ant-1]), numpy.angle( gain[n][ant-1], deg=True ) )
        str = str + " %9.0f" % s
        foutCalc.write( "%s\n" % str )

      # write 2nd line with tsys values --- #
        str = "         "
        for ant in antList :
          str = str + " %8.0f      " % tsys[n][ant-1] 
        foutCalc.write( "%s\n" % str )

      # if SEFD was nonzero, include it into the record avg --- #
        if (s > 0) :
          avg += 1./s
          navg += 1

    avgSEFD = 0.
    if (navg > 0) :
      avgSEFD = 1./(avg/navg)
    foutCalc.write("# AVG =  %9.0f\n" % avgSEFD )
    sList.append( avgSEFD )
    return sList

# avgSEFD computes SEFD averaged over a scan, for a single correlator window --- #
# ... t1 and t2 are decimal hrs since UT0, src is source name
# ... t is (npts x 1) time array, decimal hrs since UT0
# ... source is (npts x 1) list of source names
# ... gain is (npts x 15) complex array
# ... tsys is (npts x 15) real array
# ... foutCalc is file containing details of calculation (for skeptics)
# calling routine should write out info about scan number, scan UT range, correlator window to
#   file foutCalc; this routine will write out two lines (time + gains + sefd on line 1, tsys
#   on line 2) to foutCalc
#
def avgSEFD( foutCalc, t1, t2, src, t, source, gain, tsys, antList ) : 
  foutCalc.write("\n# t1,t2,src,antList = %.3f, %.3f, %s, %s\n" % (t1,t2,src,antList))

  # ... write header for foutCalc section
  foutCalc.write(" UThrs  ")
  for ant in antList :
    foutCalc.write( "  --- C%2d ---" % ant ) 
  foutCalc.write("\n")
  navg = 0
  avg = 0.

  # ... process each integration within the time interval ... #
  for n in range( 0, len(t) ) :
    if ( (t[n] >= t1) and (t[n] <= t2) and (src == source[n]) ) :
      s = SEFD( gain[n], tsys[n], antList, foutCalc )
        # ... note that gain[n] and tsys[n] are (1 x 15) arrays

      # --- write 1 line with gains and scan SEFD --- #
      str = "%8.4f " % t[n]
      for ant in antList :
        str = str + "%8.3f %5.0f" % ( abs(gain[n][ant-1]), numpy.angle( gain[n][ant-1], deg=True ) )
      str = str + " %9.0f" % s
      foutCalc.write( "%s\n" % str )

      # --- write 2nd line with tsys values --- #
      str = "         "
      for ant in antList :
        str = str + "%8.0f      " % tsys[n][ant-1] 
      foutCalc.write( "%s\n" % str )

      # --- if SEFD was nonzero, include it into the record avg --- #
      if (s > 0) :
        avg += 1./s
        navg += 1
  avgSEFD = 0.
  if (navg > 0) :
    avgSEFD = 1./(avg/navg)
  foutCalc.write("# AVG = %9.0f\n" % avgSEFD )
  return avgSEFD


# --- print phase vs time for each ant --- #
def phase( time, gain, outfile ) :
  fout = open( outfile, "w" )
  for n in range(0, len(time) ) :
    angles = numpy.array_str( numpy.angle(gain[n], deg=True ), precision=0, suppress_small=True, max_line_width=200 )
    fout.write("%8.0f %s\n" % (3600.*dechrs(time[n]), angles.strip('[,]') ) )

# --- convert hr:min:sec to decimal hrs --- #
def dechrs( hhmmss ) :
  hr,min,sec = hhmmss.split(":")
  return float(hr) + float(min)/60. + float(sec)/3600.

# --- dump out gains vs elev to look for focus problems --- #
def gainVsElev( elev, gain, outfile ) :
  fout = open( outfile, "w" )
  for n in range( 0, len(elev) ) :
    for m in range(0,15) :
      fout.write("%7.2f %5.3f" % (elev[n][m], numpy.abs(gain[n][m]) ) )
    fout.write("\n")
  fout.close()

# --- return source flux in Jy, used to predict CD or FD amps --- #
def sourceFlux( src ) :
  flux = { '3C84'     : 8.40, \
           '0854+201' : 5.00, \
           '3C273'    : 4.25, \
           '3C279'    : 16.65, \
           'M87'      : 2.20, \
           'NRAO530'  : 1.70, \
           '1921-293' : 5.15, \
           'BLLAC'    : 6.30, \
           'SGRA'     : 3.50 }
  return flux[src]
 
# --- convert day to day number --- #
def getDayNo( day ) :
  dayno = { '19mar' : '078', \
            '20mar' : '079', \
            '21mar' : '080', \
            '22mar' : '081', \
            '23mar' : '082', \
            '24mar' : '083', \
            '25mar' : '084', \
            '26mar' : '085', \
            '27mar' : '086', \
            '28mar' : '087', \
            '29mar' : '088', \
            '30mar' : '089', \
            '31mar' : '090' }
  return dayno[day]

# --- compare VLBI amps with CARMA amps --- #
# ... col = 13 for lo-band, col = 14 for hi-band
def compareAmps( VLBIfile, CARMAfile, col, fout ) :
  f1 = open( VLBIfile, "r" )
  #fout = open( "ampCmp.dat", "a" )
  fout.write("#\n")
  fout.write("# {%s col 3} vs. {%s col %2d}\n" % (VLBIfile, CARMAfile, col) )
  for line1 in f1 :
    a1 = line1.split()
    f2 = open( CARMAfile, "r" )
    for line2 in f2 :
      a2 = line2.split() 
      if ( a1[0] == a2[0] ) :  # scans match
        if ( a1[1] != a2[1] ) :  # source mismatch
          print " WARNING: scan ", a1[0], " source mismatch: ", a1[1], a2[1]
        else :                 # source match
          ampVLBI = float(a1[2])
          ampCARMA = float(a2[col-1])
          if (ampCARMA > .01) :
            fout.write(" %s  %8s  %6s %6s  %6.3f\n" % (a1[0], a1[1], a1[2], a2[col-1], ampVLBI/ampCARMA) )
          else :
            fout.write(" %s  %8s  %6s %6s    0.0\n" % (a1[0], a1[1], a1[2], a2[col-1]) )
    f2.close()
  f1.close()
  #fout.close()

# --- compute cumulative phasing efficiency for one file --- #
def cumEffic( col, infile, outfile ) :
  effic = []
  fin = open( infile, "r" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      if float(a[col-1]) > .01  :
        effic.append(a[col-1])
  fin.close()
  e = numpy.array( effic, dtype=float )
  esort = numpy.sort(e)
  print esort
  fout = open( outfile, "w" )
  fout.write("lim 0 1 0 1.05\n")
  fout.write("move 0 0\n")
  y = 0.
  n = 0
  for value in esort :
    fout.write("draw %6.3f %6.4f\n" % (value, y)  )
    n += 1
    y = float(n)/len(esort)
    fout.write("draw %6.3f %6.4f\n" % (value, y)  )
  fout.write("draw %6.3f %6.4f\n" % (1.0, y) )
  fout.close() 


# ... below version valid for 2013, not 2015; comment it out to avoid confusion when editing
# --- generate calibrations for one day (e.g., day='day075')--- #
# ... cpList is list of antennas that are phased 
# ... tsysL0,tsysR0 are the 'lo-band' systemps (from downconverters 5 and 6)
# ... tsysL1,tsysR1 are the 'hi-band' systemps (from downconverters 7 and 8)

#def oneday( day, cpList=antList ) :
#
#  # ... retrieve time and gain arrays
#  [time,gain] = getGains( day+"/wide.cal" )
#
#  # ... create decimal time array 
#  t = numpy.empty( len(time), dtype=float )
#  for n in range (0, len(time) ) :
#    t[n] = dechrs( time[n] )
#
#  # ... retrieve tsys, rmspath, tau230, elev arrays
#  tsysL0 = getVar15( day+"/tsys13L.log", t )
#  tsysR0 = getVar15( day+"/tsys13R.log", t )
#  tsysL1 = getVar15( day+"/tsys15L.log", t )
#  tsysR1 = getVar15( day+"/tsys15R.log", t )
#  rmspath = getVar( day+"/rmspath", t )
#  tau230 = getVar( day+"/tau230", t )
#  source = getSource ( day+"/uvindex.log", t )
#  elevarray = getVar15( day+"/elev", t )
#
#  # ... use elev of C1 as the elevation
#  elev = numpy.empty( (len(time)), dtype=float )
#  for n in range (0, len(time) ) :
#    elev[n] = elevarray[n][0]
#
#  # ... all array lengths should match, otherwise this routine will crash!
#  print len(time), len(gain), len(tsysL0), len(tsysR0), len(tsysL1), \
#    len(tsysR1), len(rmspath), len(tau230), len(source), len(elev)
#
#  # ... create absgain and default gain files for future use
#  absgain = numpy.empty( (len(time),15), dtype=complex )
#  defgain = numpy.empty( (len(time),15), dtype=complex )
#  jyperk = numpy.array( [1.10, 1.11, 1.00, 1.17, 1.20, 1.06, 0.93, 0.98, 1.00, 0.99, 1.25, 1.03, 0.95, 0.96, 0.98] )
#    # these are the median values from gplist on 21mar
#
#  for n in range (0, len(time)) : 
#    for i in range (0,14) :
#      absgain[n][i] = numpy.abs( gain[n][i] ) + 1j * 0.
#      defgain[n][i] = jyperk[i] * (1. + 0j)
#  print "\ndefgains:\n",defgain
#
#  # ... make table of gain vs elev for plotting
#  gainVsElev( elevarray, gain, "gainVsElev" )
#
#  # ... process the 'summ' file which lists the schedule
#  fin = open( day+"/summ", "r" )
#  fout = open( "results", "w" )
#  fout.write("#   scan        src      UTstart    UTstop   el  tau  path   CL-lo    CR-lo    FL-lo    FL-hi    FR-lo    FR-hi   pheff\n")
#  fout.write("#                                           deg        um   SEFD-Jy  SEFD-Jy  SEFD-Jy  SEFD-Jy  SEFD-Jy  SEFD-Jy\n")
#  fout.write("#\n")
#  for line in fin :
#
#    # ... print scanname, src, UTstart, UTstop (from summ file)
#    if not line.startswith("#") :
#      a = line.split()
#      src = a[0]
#      utStart = a[1]
#      utStop = a[2]
#      
#      scanname = a[7]
#      print " "
#      print scanname, src
#      fout.write(" %9s  %9s   %8s  %8s  " % (scanname,src,utStart,utStop) ) 
#
#      # ... create decimal start and stop times, print tau230, rmspath
#      t1 = dechrs( utStart ) - .0028	 # start - 10 sec
#      t2 = dechrs( utStop ) + .0028      # stop + 10 sec
#      elevavg = avgVar( t1, t2, t, elev )
#      tau230avg = avgVar( t1, t2, t, tau230 )
#      rmspathavg = avgVar( t1, t2, t, rmspath )
#      fout.write(" %2d %5.2f %4.0f" % (elevavg,tau230avg,rmspathavg ) ) 
#
#
#      # ... print SEFD for C1 L and R lo band only (hi not recorded)
#      sL0c1 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, defgain, defgain, tsysL0, [1] ) 
#      sR0c2 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, defgain, defgain, tsysR0, [1] ) 
#
#      # ... print SEFDs for phased array L and phased array R, lo and hi bands
#      sL0 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysL0, cpList ) 
#      sL0_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysL0, cpList ) 
#
#      sL1 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysL1, cpList ) 
#      #sL1_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysL1, cpList ) 
#
#      sR0 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysR0, cpList ) 
#      #sR0_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysR0, cpList ) 
#
#      sR1 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysR1, cpList ) 
#      #sR1_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysR1, cpList ) 
#
#    # all efficiencies are the same, so only write out the L-lo one
#      fout.write("   %4.2f "  % efficiency( sL0, sL0_ideal ) )
#
#      #... print predicted lo- and hi-band CARMA-CARMA correlations
#      #if (s10*sL0) > 0.1 :
#      #  amplo = 1.e4 * sourceFlux(src) / math.sqrt(s10*sL0)
#      #else :
#      #  amplo = 0.
#      #fout.write("  %6.3f" % amplo )
#
#      fout.write("\n")
#
#  # ... finished
#  fin.close()
#  fout.close()


# --- compute improvement, relative to one 10-m telescope, for n10 and n6 phased ants --- #
def idealSEFD(n10, n6) :
  tsys = 100.
  g10 = 65.
  g6 = 145.
  merit1 = 1./(g10 * tsys)
  merit2 = 1./( (n10 + n6) * tsys ) * pow(n10/math.sqrt(g10) + n6/math.sqrt(g6), 2.)
  merit3 = 1./tsys * pow(n10/g10 + n6/g6, 2.)/(n10/g10 + n6/g6)
  print "equal power weighting: %5.3f" %  (merit2/merit1)
  print "optimum weighting: %5.3f" %  (merit3/merit1)


# reads phs.hisotry file
def readHist( infile, outfile ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  for line in fin:
    if ( not line.startswith("#") ) :
      a = line.split()
      for n in range(0,15) :
        phs[n] = math.pi/180. * float(a[n+6])
  fin.close()


# calculation of phasing effic for Heino Falcke
# procedure: (1) fit passband, rewrite to apply; (2) selfcal with interval=10000 to
#   remove any average phase offsets (these would be taken out in the beamformer);
#   rewrite data once again to apply these; (3) now do record-by-record selfcal;
#   we will assume that these selfcal solutions would have given perfect phasing
#   (not true in cases of poor s/n, but seems roughly OK for SgrA); since we were
#   running the phasing script during these observations, these gains are the
#   residual errors; phasing effic = vector amp of these gains (all with amp=1,
#   all assumed to have same tsys), vs scalar average
# phasing effic = vectorSum/scalarSum, but scalarSum is just the number of antennas

def basicEffic( antList, infile, outfile ) :
  [time,gain] = getGains( infile )
  fout = open(outfile, "w")
  for n in range (0, len(time) ) :
    vectorSum = 0.+0.j
    scalarSum = 0.
    for ant in antList :
      vectorSum = vectorSum + gain[n,ant-1]
      scalarSum = scalarSum + 1. 
    fout.write("%s  %7.4f  %6.4f\n" % (time[n], dechrs(time[n]), numpy.abs(vectorSum)/scalarSum ))
  fout.close()

def old( ) :
  basicEffic( antList, "day075/s.vis", "day075/s.phEffic" )
  cumEffic( 3, "day075/s.phEffic", "day075cum.wip" ) 
  basicEffic( antList, "day080/s.vis", "day080/s.phEffic" )
  cumEffic( 3, "day080/s.phEffic", "day080cum.wip" ) 
  basicEffic( antList, "day081/s.vis", "day081/s.phEffic" )
  cumEffic( 3, "day081/s.phEffic", "day081cum.wip" ) 

# .... valid for 2013; comment out to avoid confusion
# --- generate record-by-record gains for one day (e.g., day='day075')--- #
# ... cpList is list of antennas that are phased 
# ... tsysL0,tsysR0 are the 'lo-band' systemps (from downconverters 5 and 6)
# ... tsysL1,tsysR1 are the 'hi-band' systemps (from downconverters 7 and 8)

#def onedayRbyR( day, cpList=antList ) :
#
#  # ... retrieve time and gain arrays
#  [time,gain] = getGains( day+"/wide.cal" )
#
#  # ... create decimal time array 
#  t = numpy.empty( len(time), dtype=float )
#  for n in range (0, len(time) ) :
#    t[n] = dechrs( time[n] )
#
#  # ... retrieve tsys, rmspath, tau230, elev arrays
#  tsysL0 = getVar15( day+"/tsys13L.log", t )
#  tsysR0 = getVar15( day+"/tsys13R.log", t )
#  tsysL1 = getVar15( day+"/tsys15L.log", t )
#  tsysR1 = getVar15( day+"/tsys15R.log", t )
#  rmspath = getVar( day+"/rmspath", t )
#  tau230 = getVar( day+"/tau230", t )
#  source = getSource ( day+"/uvindex.log", t )
#  elevarray = getVar15( day+"/elev", t )
#
#  # ... use elev of C1 as the elevation
#  elev = numpy.empty( (len(time)), dtype=float )
#  for n in range (0, len(time) ) :
#    elev[n] = elevarray[n][0]
#
#  # ... all array lengths should match, otherwise this routine will crash!
#  print len(time), len(gain), len(tsysL0), len(tsysR0), len(tsysL1), \
#    len(tsysR1), len(rmspath), len(tau230), len(source), len(elev)
#
#  # ... create absgain and default gain files for future use
#  absgain = numpy.empty( (len(time),15), dtype=complex )
#  defgain = numpy.empty( (len(time),15), dtype=complex )
#  jyperk = numpy.array( [1.10, 1.11, 1.00, 1.17, 1.20, 1.06, 0.93, 0.98, 1.00, 0.99, 1.25, 1.03, 0.95, 0.96, 0.98] )
#    # these are the median values from gplist on 21mar
#
#  for n in range (0, len(time)) : 
#    for i in range (0,14) :
#      absgain[n][i] = numpy.abs( gain[n][i] ) + 1j * 0.
#      defgain[n][i] = jyperk[i] * (1. + 0j)
#  print "\ndefgains:\n",defgain
#
#  # ... make table of gain vs elev for plotting
#  gainVsElev( elevarray, gain, "gainVsElev" )
#
#  # ... process the 'summ' file which lists the schedule
#  fout = open( "RbyR", "a" )
#  fin = open( day+"/summ", "r" )
#
#  for line in fin :
#    if not line.startswith("#") :
#      a = line.split()
#      src = a[0]
#      utStart = a[1]
#      utStop = a[2]
#      scanname = a[7]
#      print " "
#      print scanname, src
#      fout.write("#\n")
#      fout.write("#    day %s, scan %s, source %s, UT %8s -- %8s \n" % \
#           ( getDayNo(day), scanname, src, utStart, utStop) ) 
#      fout.write("#\n#    UT        UThrs     el  tau  path      CL-lo    CR-lo    FL-lo    FL-hi    FR-lo    FR-hi   pheff\n")
#
#      for n in range(0,len(time)) :
#
#        if (t[n] > dechrs( utStart ) - .0028) and (t[n] < dechrs( utStop ) + .0028) :
#          fout.write("  %8s %10.6f    %2d %5.2f %4.0f   " % (time[n],t[n],elev[n],tau230[n],rmspath[n] ) ) 
#
#          sL0c1 = SEFD( defgain[n], tsysL0[n], [1] )
#          fout.write("%8.0f " % (100.*round(sL0c1/100.)) )
#
#          sR0c1 = SEFD( defgain[n], tsysR0[n], [1] )
#          fout.write("%8.0f " % (100.*round(sR0c1/100.)) )
#
#          sL0 = SEFD( gain[n], tsysL0[n], cpList )
#          fout.write("%8.0f " % (100.*round(sL0/100.)) )
#
#          sR0 = SEFD( gain[n], tsysR0[n], cpList )
#          fout.write("%8.0f " % (100.*round(sR0/100.)) )
#
#          sL1 = SEFD( gain[n], tsysL1[n], cpList )
#          fout.write("%8.0f " % (100.*round(sL1/100.)) )
#
#          sR1= SEFD( gain[n], tsysR1[n], cpList )
#          fout.write("%8.0f " % (100.*round(sR1/100.)) )
#
#          sL0_ideal = SEFD( absgain[n], tsysL0[n], cpList ) 
#          fout.write("   %4.2f\n"  % efficiency( sL0, sL0_ideal ) )
#
#      # now print out average values, as before
#      t1 = dechrs( utStart ) - .0028	 # start - 10 sec
#      t2 = dechrs( utStop ) + .0028      # stop + 10 sec
#      elevavg = avgVar( t1, t2, t, elev )
#      tau230avg = avgVar( t1, t2, t, tau230 )
#      rmspathavg = avgVar( t1, t2, t, rmspath )
#      fout.write("# ------- AVG -------   ")
#      fout.write(" %2d %5.2f %4.0f   " % (elevavg,tau230avg,rmspathavg ) ) 
#
#      # ... print SEFD for C1 L and R lo band only (hi not recorded)
#      sL0c1 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, defgain, defgain, tsysL0, [1] ) 
#      sR0c2 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, defgain, defgain, tsysR0, [1] ) 
#
#      # ... print SEFDs for phased array L and phased array R, lo and hi bands
#      sL0 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysL0, cpList ) 
#      sL0_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysL0, cpList ) 
#
#      sL1 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysL1, cpList ) 
#      #sL1_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysL1, cpList ) 
#
#      sR0 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysR0, cpList ) 
#      #sR0_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysR0, cpList ) 
#
#      sR1 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysR1, cpList ) 
#      #sR1_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysR1, cpList ) 
#
#    # all efficiencies are the same, so only write out the L-lo one
#      fout.write("   %4.2f "  % efficiency( sL0, sL0_ideal ) )
#      fout.write("\n")
#
#  fin.close()  
###  fout.close()

def makeScanList( day ) :
  scanname = []
  utStart = []
  utStop = []
  fin = open( day+"/summ", "r" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      utStart.append( a[1l])
      utStop.append( a[2] )
      scanname.append( a[7] )
  fin.close()
  return [scanname, utStart, utStop ]
  
def findscan( t, scanname, utStart, utStop ) :
  sname = " "
  ut1 = " "
  ut2 = " "
  for n in range(0, len(utStart) ) :
    print " n = %d" % n
    if (t > dechrs( utStart[n] ) - .0028) and (t < dechrs( utStop[n] ) + .0028) :
      sname = scanname[n]
      ut1 = utStart[n]
      ut2 = utStop[n]
      break
  return [ sname, ut1, ut2 ]
  

#def doit2013() :
#   
#  antList = [2,3,4,5,6,8,9,14]    # use for 21mar and 22mar
#  #onedayRbyR( "21mar", cpList=antList ) 
#  #onedayRbyR( "22mar", cpList=antList ) 
#  antList = [2,4,5,6,8,9,13,14]   # use for 23mar and beyond
#  #onedayRbyR( "23mar", cpList=antList ) 
#  #onedayRbyR( "25mar", cpList=antList ) 
#  #onedayRbyR( "26mar", cpList=antList ) 
#  onedayRbyR( "27mar", cpList=antList ) 
  
# ========== added 24feb2015 ========= #

# run from subdirectory of vlbi/2011
def flux2011( ) :
  time = []
  flux = []
  fin = open( "fluxes.txt", "r" )
  fout = open( "fluxVsTime", "w" )

  for line in fin :
    a = line.split()
    if len(a) == 2 :
      time.append( 24.*float(a[0]) )  # dechrs
      flux.append( float(a[1]) )      # flux

  # ... create decimal time array 
  t = numpy.array( time, dtype=float )
  print t

  # ... retrieve tsys, rmspath, tau230, elev arrays
  rmspath = getVar( "rmspath.log", t )
  tau230 = getVar( "tau230.log", t )
  source = getSource ( "uvindex.log", t )
  elevarray = getVar15( "elev.log", t )

  # ... use elev of C1 as the elevation
  elev = numpy.empty( (len(time)), dtype=float )
  for n in range (0, len(time) ) :
    elev[n] = elevarray[n][0]

  # ... all array lengths should match, otherwise this routine will crash!
  print len(time), len(flux), len(rmspath), len(tau230), len(source), len(elev)

  fout.write("#   UThrs     S(Jy)    el    tau   path   source\n")
  for n in range (0, len(time)) : 
    fout.write("%10.5f   %6.3f   %4.1f   %4.2f  %4.0f   %s\n" % (t[n], flux[n], elev[n], tau230[n], rmspath[n], source[n] ))
  fin.close()
  fout.close()

# try to understand what time points are missing from flux list
def diagnose() :
  fin = open( "fluxes.txt", "r" )
  fout = open( "junk", "w" )
  for line in fin :
    a = line.split()
    if len(a) == 2 :
      time = 24.*float(a[0])   # dechrs
      h = int(time) 
      m = int( 60.*(time-float(h) ) )
      s = int( 60.*(60.*(time-float(h)) - float(m)))
      print h, m, s
      fout.write("%3d:%2d:%2d\n" % (h,m,s) )
  fin.close()
  fout.close()

# this routine was written specifically to analyze fluxes from flux2.csv, 2015 data
def fluxSummary() :
  fin = open("flux2.csv", "r" )
  goodcols = [9,10,11,12,15,17,18,20,22,23,25]
	# use fluxes in these columns from spreadsheet
  for line in fin :
    a = line.split(",")
    src = a[0]
    S = []
    sum = 0.
    n = 0
    for j in goodcols :
      if (len(a[j]) == 0) or (a[j] == " ") :
        S.append(  0. )
      else :
        S.append( float(a[j]) )
        sum = sum + float(a[j])
        n = n+1
    avg = 0.
    if n > 0 :
      avg = sum/n
    print "%10s  %.2f  %d" % (src, avg, n)
    #if n > 2 :
    #  scale[i] = numpy.array(S)/avg
    #  print scale
  fin.close()
   
# this routine written to diagnose effect of wild gain fluctuations from C14

def c14( ) :
    fout = open( "c14L.log", "w" )
    time, gainComplex = getGains( "wide.av" )

  # create decimal time array 
    t = numpy.empty( len(time), dtype=float )
    for n in range (0, len(time) ) :
      t[n] = dechrs( time[n] )

    tsys = getVar15( "tsysL.log", t ) 
    for n in range(0, len(time) ) :
      g = abs(gainComplex[n][13])
      ts = tsys[n][13]
      fout.write("%s  %.5f  %.0f  %.3f  %.3f  %.3f  %.3f  %.0f\n" % \
        (time[n],t[n],ts, g, g*g, ts*g, ts*g*g, SEFD( gainComplex[n], tsys[n], [14] ) ) ) 
    fout.close()


# --- generate record-by-record SEFDs and scan-averaged SEFDs for one day in 2015 experiment
# ... cpList is list of antennas that are phased (always the same in 2015)
# ... day is a date, like 31mar
# ... works with: gainL = selfcal gains for LCP, uvrange(20,1000), all sources including SgrA
# ...             gainR = selfcal gains for RCP, uvrange(20,1000), all sources including SgrA
#                 TsysL5,tsysR5, etc
# t1str,t2str is time range to use to derive SgrA gains
#
def oneday2015( day, t1str, t2str, C1gain=None, cpList=[2,3,4,5,6,13,14,15], makeGainPlots=False ) :

  # retrieve time and gain arrays
    print "... retrieving L gains"
    [timeL,gainLtmp] = getGains( "wideL.vlb" )
    print "... retrieving R gains"
    [timeR,gainRtmp] = getGains( "wideR.vlb" )

  # drop any time,gain pairs where both L and R gains are not present
    print "... number of L gains = %d; number of R gains = %d" % ( len(timeL),len(timeR) )
    time = []
    iL = 0
    iR = 0
    Lgain = []
    Rgain = []
    while (iL < len(timeL)) and (iR < len(timeR)) :
      if timeL[iL] == timeR[iR] :
        time.append( timeL[iL] )
        Lgain.append( gainLtmp[iL] )
        Rgain.append( gainRtmp[iR] )
        iL = iL + 1
        iR = iR + 1
      else :
        tL = dechrs( timeL[iL] )
        tR = dechrs( timeR[iR] )
        if tL > tR :
          print "...  discarding Rgain for time %s" % timeR[iR]   
          iR = iR + 1
        else :
          print "...  discarding Lgain for time %s" % timeL[iL]   
          iL = iL + 1
    gainLtmp = numpy.array( Lgain )
    gainRtmp = numpy.array( Rgain )
    print "... after reconciliation, number of L gains = %d; number of R gains = %d" % ( len(gainLtmp),len(gainRtmp) )
    
  # open output files, label with current date,time so we can keep track of old vs new ones
    foutCalc = open( "SEFD_calc."+day, "a" )   # 2 lines per integration, one band; gain, Tsys, SEFD for each record
    foutRbyR = open( "SEFD_RbyR."+day, "a" )	 # one line per integration, all bands; SEFD values only
    foutAvg = open( "SEFD_avg."+day, "a" )     # one line per scan, all bands; scan-averaged values only
    foutSummary = open( "SEFD_summary."+day, "a" )     # one line per scan, Cr and Cp avgSEFD; Cp 10-90 percentile SEFDs

    nowStr = datetime.datetime.now().strftime("%m/%d/%Y %H:%M")
    foutAvg.write("# [vCal output %s]\n#\n" % nowStr )
    foutCalc.write("# [vCal output %s]\n" % nowStr )
    foutRbyR.write("# [vCal output %s]\n" % nowStr )
    foutSummary.write("# [vCal output %s]\n" % nowStr )
    foutAvg.write("# day scan    source   UTstart  el   tau path      Cr-L     Cr-R     Cp-1L    Cp-1R    Cp-3L    Cp-3R    Cp-5L    Cp-5R    Cp-7L    Cp-7R   pheff\n")

  # rescale gain amps on SgrA, fill in gains on C1 during scans (note: SgrGain reads source.log and axisrms.log)
    print "... fill in SgrA and C1 Lgains"
    plotName = None
    if makeGainPlots :
      plotName = "Lgains.%s.pdf" % day
    foutCalc.write("#\n# Lgains\n" )
    gainL = SgrGain( time, gainLtmp, t1str, t2str, C1gain, plotName, foutCalc )
    print "... fill in SgrA and C1 Rgains"
    if makeGainPlots :
      plotName = "Rgains.%s.pdf" % day
    foutCalc.write("#\n# Rgains\n" )
    gainR = SgrGain( time, gainRtmp, t1str, t2str, C1gain, plotName, foutCalc )

  # create decimal time array 
    t = numpy.empty( len(time), dtype=float )
    for n in range (0, len(time) ) :
      t[n] = dechrs( time[n] )

  # retrieve tsys, rmspath, tau230, elev arrays; note that only values corresponding exactly to "t" array are returned
    print "... retrieve system temperatures"
    tsysL5 = getVar15( "tsys5L.log", t )
    tsysR5 = getVar15( "tsys5R.log", t )
    tsysL6 = getVar15( "tsys6L.log", t )
    tsysR6 = getVar15( "tsys6R.log", t )
    tsysL7 = getVar15( "tsys7L.log", t )
    tsysR7 = getVar15( "tsys7R.log", t )
    tsysL8 = getVar15( "tsys8L.log", t )
    tsysR8 = getVar15( "tsys8R.log", t )
    print "... retrieve rmspath, tau230, sourceName, elevation"
    rmspath = getVar( "rmspath.log", t )
    tau230 = getVar( "tau230.log", t )
    source = getSource ( "source.log", t )
    elevarray = getVar15( "elev.log", t )
    elev = elevarray[:,7]
      # assume C8 elev applies to the whole array

  # all array lengths should match, otherwise this routine will crash!
    print "... array lengths: ", len(time), len(gainL), len(gainR), len(tsysL5), len(tsysR5), len(tsysL6), \
      len(tsysR6), len(rmspath), len(tau230), len(source), len(elev)

  # create absgain array for future use in computing phasing efficiency
    absgain = numpy.abs( gainL )

  # make table of gain vs elev for plotting
    gainVsElev( elevarray, gainL, "gainVsElev" )

  # process the 'summ' file which lists the schedule
    print "... process summ file"
    fin = open( "summ", "r" )


  # process each scan
    for line in fin :
      if not line.startswith("#") :
        a = line.split()
        if a[4] != "vlb" :
          print "ERROR - col 5 is not 'vlb'"
        src = a[0]
        utStart = a[1]
        utStop = a[2]
        scanname = a[7]

        print " "
        print scanname, src

        foutCalc.write("#\n# =============================================================================================================================\n")
        foutCalc.write("# day %s, scan %s, source %s, UT %8s -- %8s \n" % \
           ( getDayNo(day), scanname, src, utStart, utStop) ) 

        foutRbyR.write("#\n")
        foutRbyR.write("#    day %s, scan %s, source %s, UT %8s -- %8s \n" % \
         ( getDayNo(day), scanname, src, utStart, utStop) ) 
        foutRbyR.write("#    UT      UThrs     el   tau path       Cr-L     Cr-R    Cp-1L    Cp-1R    Cp-3L    Cp-3R    Cp-5L    Cp-5R    Cp-7L    Cp-7R   pheff\n")

        t1 = dechrs( utStart ) - .0028	 # start - 10 sec
        t2 = dechrs( utStop ) + .0028      # stop + 10 sec
   
      # generate [time x band] array of record-by-record SEFDs for each recorder; last entry is scan average
        ARbyR = []
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cr-L', t, gainL, tsysL5, [1] ) )
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cr-R', t, gainR, tsysL5, [1] ) )
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-1L', t, gainL, tsysL5, cpList ) )
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-1R', t, gainR, tsysR5, cpList ) )
          # comment out for 21,22mar2015
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-3L', t, gainL, tsysL6, cpList ) )
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-3R', t, gainR, tsysR6, cpList ) )
          # comment out for 21,22mar2015
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-5L', t, gainL, tsysL7, cpList ) )
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-5R', t, gainR, tsysR7, cpList ) )
          # comment out for 21,22mar2015
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-7L', t, gainL, tsysL8, cpList ) )
        ARbyR.append( AscanSEFD( foutCalc, t1, t2, 'Cp-7R', t, gainR, tsysR8, cpList ) )
          # comment out for 21,22mar2015
        
      # compute phasing efficiency for band6L, assume it is about the same for other bands
        idealSEFD = AscanSEFD( foutCalc, t1, t2, 'Cp-3L-ideal', t, absgain, tsysL6, cpList )

      # write record-by-record summary, all bands on one line     
        i = 0    # index for ARbyR arrays, which cover only (t1,t2), hence are shorter
        for n in range(0,len(time)) :
          if (t[n] > t1) and (t[n] < t2 ) :     # valid integration in this scan
            foutRbyR.write("  %8s %8.4f    %2d %5.2f %4.0f   " % (time[n],t[n],elev[n],tau230[n],rmspath[n] ) ) 

          # this is the normal code
            for j in range(0,10) :
              foutRbyR.write("%8.0f " % (100.*round(ARbyR[j][i]/100.)) )

          # code below is ony for 21,22mar2015
            # for j in range(0,2) :
            #   foutRbyR.write("%8.0f " % (100.*round(ARbyR[j][i]/100.)) )
            # for j in range(2,6) :
            #   foutRbyR.write("%8.0f " % (100.*round(ARbyR[j][i]/100.)) )
            #   foutRbyR.write("    ---- " )

            foutRbyR.write("   %4.2f\n"  % efficiency( ARbyR[4][i], idealSEFD[i] ) )
            i = i + 1

      # finally, print line with average values into foutRbyR and foutAvg
        elevavg = avgVar( t1, t2, t, elev )
        tau230avg = avgVar( t1, t2, t, tau230 )
        rmspathavg = avgVar( t1, t2, t, rmspath )
        foutRbyR.write("# ------ AVG ------   ")
        foutRbyR.write(" %2d %5.2f %4.0f   " % (elevavg,tau230avg,rmspathavg ) ) 
        foutAvg.write("  %3s %4s  %8s  %8s" % ( getDayNo(day), scanname, src, utStart) )
        foutAvg.write("  %2d %5.2f %4.0f   " % (elevavg,tau230avg,rmspathavg ) )  
        foutSummary.write("  %3s %4s  %8s  %8s" % ( getDayNo(day), scanname, src, utStart) )
        foutSummary.write("  %2d %5.2f %4.0f   " % (elevavg,tau230avg,rmspathavg ) )  

      # this is the normal code
        for j in range(0,10) :
          foutRbyR.write("%8.0f " % (100.*round(ARbyR[j][-1]/100.)) )
          foutAvg.write("%8.0f " % (100.*round(ARbyR[j][-1]/100.)) )

      # this is the code for 21,22 mar 2015 when RCP beamformers were misphased  
        #for j in range(0,2) :
        #  foutRbyR.write("%8.0f " % (100.*round(ARbyR[j][-1]/100.)) )
        #  foutAvg.write("%8.0f " % (100.*round(ARbyR[j][-1]/100.)) )
        #for j in range(2,6) :
        #  foutRbyR.write("%8.0f " % (100.*round(ARbyR[j][-1]/100.)) )
        #  foutRbyR.write("    ---- " )
        #  foutAvg.write("%8.0f " % (100.*round(ARbyR[j][-1]/100.)) )
        #  foutAvg.write("    ---- " )

        foutRbyR.write("   %4.2f\n"  % efficiency( ARbyR[4][-1], idealSEFD[-1] ) )
        foutAvg.write("   %4.2f\n"  % efficiency( ARbyR[4][-1], idealSEFD[-1] ) )

      # for plot of SEFD range, reorder by SEFD, then find 10 and 90 percentile values
        meanCr = 0.5 * ( ARbyR[0][-1] + ARbyR[1][-1] )
          # avg Cr SEFD for L and R, over the entire scan - one number
        meanSEFD = numpy.mean( ARbyR[2:], axis=0 ) 
          # avg Cp SEFD for the 8 windows (4 L, 4 R) for each 10 sec record - an array
        sortedCpSEFD = numpy.sort( meanSEFD[0:-1] )
          # meanSEFD sorted from lowest to highest, omitting last value, which is the avg
        if len(sortedCpSEFD) > 2 :
          nn1 = int( round(0.1 * len(sortedCpSEFD)) )
          nn2 = int( round(0.9 * len(sortedCpSEFD)) )
          print "%d %d  %.0f  %.0f %.0f %.0f" % ( nn1,nn2,meanCr,meanSEFD[-1],sortedCpSEFD[nn1-1], sortedCpSEFD[nn2-1])
          foutSummary.write(" %8.0f  %8.0f %8.0f %8.0f    %6.4f   %6.4f %6.4f %6.4f\n" % \
            (meanCr, meanSEFD[-1], sortedCpSEFD[nn1-1], sortedCpSEFD[nn2-1], numpy.log10(meanCr), \
            numpy.log10(meanSEFD[-1]), numpy.log10(sortedCpSEFD[nn1-1]), numpy.log10(sortedCpSEFD[nn2-1]) ) )
              # printing logs for possible logarithmic plot
        else :
          foutSummary.write(" 0. 0. 0. 0. 0. 0. 0. 0.\n")


    fin.close()  
    foutCalc.close()
    foutRbyR.close()
    foutAvg.close()
    foutSummary.close()

def makeSEFDplot( infile, outfile ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  for line in fin :
    a = line.split()
    col = 1
    if (a[2] == "SGRA" ) : col = 2
    if (a[2] == "M87" ) : col = 4
    fout.write("move %s %s\n" % (a[11],a[1]) )
    fout.write("color %d\n" % col)
    fout.write("symbol 0\n")
    fout.write("dot\n")
    fout.write("move %s %s\n" % (a[12],a[1]) )
    fout.write("symbol 16\n" )
    fout.write("dot\n")
    fout.write("move %s %s\n" % (a[13],a[1]) )
    fout.write("draw %s %s\n" % (a[14],a[1]) )
    fout.write("expand .8\n")
    fout.write("move 6.1 %.2f\n" % (float(a[1])+.5) )
    fout.write("putlabel 0 %s\n" % a[1])
    fout.write("move 6.5 %.2f\n" % (float(a[1])+.5) )
    fout.write("putlabel 0 %s\n" % a[3])
    fout.write("move 7.2 %.2f\n" % (float(a[1])+.5) )
    fout.write("putlabel 0 %s\n" % a[2])
    fout.write("expand 1\n")
    fout.write("color 1\n")
  fin.close()
  fout.close()
   
