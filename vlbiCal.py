# vlbiCal.py
#
# computes SEFD for phased and unphased antennas during vlbi run

import math
import time
import numpy
import subprocess
import shlex

antList = [2,3,4,5,8,12,13,14]
 
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
      gainComplex[ng][nant-1] = float(a[5]) + 1j * float(a[6])
    elif ( len(a) > 0 ) :
      nant = int(a[1]) 
      if ( nant < 16 ) :
        gainComplex[ng][nant-1] = float(a[4]) + 1j * float(a[5])
  return [time, gainComplex]
   
# --- this is a helper file that parses output of gplist --- #
def parseGains( input ) :
  time = []
  lines = input.split("\n")
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
      gainComplex[ng][nant-1] = float(a[5]) + 1j * float(a[6])
    elif ( len(a) > 0 ) :
      nant = int(a[1]) 
      if ( nant < 16 ) :
        gainComplex[ng][nant-1] = float(a[4]) + 1j * float(a[5])
  return [time, gainComplex]


# --- compute R/L gain ratio from gplist outputs Rgains and Lgains --- #
# ... selfcal vis=wide.cal select=pol(RR) options=amplitude,apriori,noscale ...
# ... gplist vis=wide.cal options=all > Rgains
# ... selfcal vis=wide.cal select=pol(LL) options=amplitude,apriori,noscale ...
# ... gplist vis=wide.cal options=all > Lgains

def RLratio( RgainsFile, LgainsFile, outfile ) :
  ratios = []
  fin = open( RgainsFile, "r" )
  [timeR, gainComplexR] = parseGains( fin.read() )
  fin.close()
  fin = open( LgainsFile, "r" )
  [timeL, gainComplexL] = parseGains( fin.read() )
  fin.close()
  fout = open( outfile, "w" )
  if len(timeR) != len(timeL) :
    print "len(timeR) = %d != len(timeL) = % d - exiting" % ( len(timeR),len(timeL) )
    return
  for n in range( 0, len(timeR) ) :
    if (timeR[n] != timeL[n] ) :
      print timeR[n], " != ", timeL[n], " - exiting"
      return
    else :
      ratio = (abs(gainComplexL[n])*abs(gainComplexL[n])) / (abs(gainComplexR[n])*abs(gainComplexR[n]) )
      fout.write(" %s %s\n" % (timeR[n],numpy.array_str( ratio, precision=3, max_line_width=200 ) ) )
      ratios.append( ratio )
  maskedRatios = numpy.ma.array( ratios, mask=(numpy.isnan(ratios) ) )	 # mask off the nans
  print "AVG: ", numpy.array_str( numpy.ma.mean(maskedRatios, axis=0 ), precision=3, max_line_width=200 )
  fout.close()
         
  
  
# --- read varplt log containing tsys or elev, 6 lines per vis record, 4 values per line --- #
# ... t is the decimal time array
def getVar15( infile, t ) :
  varlist = []
  fin = open( infile, "r" )
  n1 = 0
  nline = 0
  for line in fin:
    if ( not line.startswith("#") ) :
      a = line.split()
      if ( len(a) == 6) : 
        print a[1], t[nline]
        if ( (nline > len(t)) or (abs(dechrs(a[1]) - dechrs(t[nline]) ) > .0003)) :
          print "skipping %s record at UT %s - time not in gains array" % (infile,a[1])
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
  return numpy.reshape( vl, (-1,15) )

# --- read varplt log containing single variable like rmspath or tau225 --- "
def getVar( infile, t ) :
  var = []
  fin = open( infile, "r" )
  nline = 0
  for line in fin:
    if ( not line.startswith("#") ) :
      a = line.split()
      if ( abs(dechrs(a[1]) - t[nline]) > .0003) :
        print "skipping %s record at UT %s - time not in gains array" % (infile,a[1])
      else :
        var.append( a[2] )
        nline = nline + 1
  fin.close()
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
        if ( abs(dtime - t[nline]) > .0003) :
          print "skipping %s record at UT %s - time not in gains array" % (infile,a[1])
        else :
          source.append( a[1] )
          nline = nline + 1

# --- generate calibrations for one day (e.g., day='01apr')--- #
# ... cpList is list of antennas that are phased (used if summ file says 'CP')
def oneday( day, cpList ) :

  # ... retrieve time and gain arrays
  [time,gain] = getGains( infile )

  # ... create decimal time array 
  t = numpy.empty( len(time), dtype=float )
  for n in range (0, len(time) ) :
    t[n] = dechrs( time[n] )

  # ... retrieve tsys, rmspath, tau230, elev arrays
  tsys0 = getVar15( day+"/tsys.win13", t )
  tsys1 = getVar15( day+"/tsys.win15", t )
  rmspath = getVar( day+"/rmspath", t )
  tau230 = getVar( day+"/tau230", t )
  source = getSource ( day+"/uvindex.log", t )
  elev = getVar15( day+"/elev", t )

  # ... all array lengths should match, otherwise this routine will crash!
  print len(time), len(gain), len(tsys0), len(tsys1), len(rmspath), len(tau230), len(source), len(elev)

  # ... create absgain and default gain files for future use
  absgain = numpy.empty( (len(time),15), dtype=complex )
  defgain = numpy.empty( (len(time),15), dtype=complex )
  for n in range (0, len(time)) : 
    for i in range (0,14) :
      absgain[n][i] = numpy.abs( gain[n][i] ) + 1j * 0.
      defgain[n][i] = 1. + 0j

  # ... make table of gain vs elev for plotting
  gainVsElev( elev, gain, day+".gainVsElev" )

  # ... process the 'summ' file which lists the schedule
  fin = open( day+"/summ", "r" )
  dayno = getDayNo( day )
  fout = open( "day"+dayno+".results", "w" )
  fout.write("#   scan        src     UTstart   UTstop    el  tau  path    C1-lo     C1-hi     C4-lo     CP-hi  ph-eff     CD     FD\n")
  fout.write("#                                          deg        um    SEFD-Jy   SEFD-Jy   SEFD-Jy   SEFD-Jy           amp    amp\n")
  for line in fin :

    # ... print scanname, src, UTstart, UTstop, elev (from summ file)
    a = line.split()
    src = a[4]
    hr,min,sec = a[8].split(":")
    scanname = dayno + "-" + str(hr)+str(min)+str(sec)
    fout.write(" %10s %9s  %8s  %8s   %2s" % (scanname,src,a[8],a[9],a[6] ) ) 

    # ... create decimal start and stop times, print tau230, rmspath
    t1 = dechrs( a[8] )  # scan start time
    t2 = dechrs( a[9] )  # scan stop time
    tau230avg = avgVar( t1, t2, t, tau230 )
    rmspathavg = avgVar( t1, t2, t, rmspath )
    fout.write(" %5.2f %4.0f" % (tau230avg,rmspathavg ) ) 

    print a[8], t1, a[9], t2, src

    # ... print SEFDs for C1 lo and hi bands; use defgain if no prescan data
    s10 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, gain, defgain, tsys0, [1] ) 
    s11 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, gain, defgain, tsys1, [1] ) 

    # ... print SEFDs for C4 or CP (phased array) lo and hi bands
    if (a[0] == "C4")  :
      s40 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, gain, defgain, tsys0, [4] ) 
    else :
      s40 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsys0, cpList ) 

    if (a[1] == "C4") :
      s41 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, gain, defgain, tsys1, [4] ) 
      fout.write("   -1  " )
    else :
      s41 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsys1, cpList ) 
      s41_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsys1, cpList ) 
      if (s41 > 0.) :
        effic = s41_ideal/s41
      else :
        effic = 0.
      fout.write("  %4.2f "  % effic )

    # ... print predicted lo- and hi-band CARMA-CARMA correlations
    if (s10*s40) > 0.1 :
      amplo = 1.e4 * sourceFlux(src) / math.sqrt(s10*s40)
    else :
      amplo = 0.
    if (s11*s41) > 0.1 :
      amphi = 1.e4 * sourceFlux(src) / math.sqrt(s11*s41)
    else :
      amphi = 0.
    fout.write("  %6.3f %6.3f\n" % (amplo,amphi) )

  # ... finished
  fin.close()
  fout.close()


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

# --- wrapper for avgSEFD that uses default gains if no pre- or post-scan data --- #
def scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsys, antList ) :
  s = avgSEFD( t1, t2, src, t, source, gain, tsys, antList ) 
  if (s > 0.) :
    fout.write(" %8.0f " % s )
  else :
    s = avgSEFD( t1, t2, src, t, source, defgain, tsys, antList )  
    fout.write(" %8.0f*" % s )
  return s

# --- compute SEFD averaged over a scan --- #
# ... t1 and t2 are decimal hrs since UT0, src is source name
# ... t is (npts x 1) time array, decimal hrs since UT0
# ... source is (npts x 1) list of source names
# ... gain is (npts x 15) complex array
# ... tsys is (npts x 15) real array
def avgSEFD( t1, t2, src, t, source, gain, tsys, antList ) : 
  print "\nt1, t2, src, antList = %.4f, %.4f, %s, %s" % (t1,t2,src,antList)
  navg = 0
  avg = 0.
  for n in range( 0, len(t) ) :
    if ( (t[n] >= t1) and (t[n] <= t2) and (src == source[n]) ) :
      s = SEFD( gain[n], tsys[n], antList )
        # ... note that gain[n] and tsys[n] are (1 x 15) arrays
      print "%8.4f %9.0f" % (t[n], s)
      if (s > 0) :
        avg += 1./s
        navg += 1
  if (navg > 0) :
    print "   AVG = %9.0f" % (1./(avg/navg))
    return 1./(avg/navg)
  else :
    print "   AVG =       0"
    return 0.


# --- compute SEFD for single or phased ants, for 1 integration --- #
# ... onegain = 1 x 15 complex array
# ... onetsys = 1 x 15 real array
def SEFD( onegain, onetsys, antList ) :
  jyperk = [65.,65.,65.,65.,65.,65.,145.34,145.34,145.34,145.34,145.34,145.34,145.34,145.34,145.34]
  nants = 0
  sum = 0.
  for ant in antList :
    magg = numpy.abs( onegain[ant-1] )
    if ( (magg > 0.1) and (magg < 10.) and (onetsys[ant-1] > 10. ) ) :
      nants += 1
      sum += 1./( onegain[ant-1] * math.sqrt( jyperk[ant-1] * onetsys[ant-1] ) )
  if (nants > 0) :
    return nants/(numpy.abs(sum) * numpy.abs(sum))  
  else :
    return 0.

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
  flux = { '0854+201' : 4.39, \
           '3C273'    : 6.18, \
           '1633+382' : 3.09, \
           '3C279'    : 13.84, \
           '1749+096' : 3.98, \
           'M87'      : 1.66, \
           'NRAO530'  : 3.32, \
           '1921-293' : 5.41, \
           '3C454.3'  : 11.42, \
           '3C345'    : 2.14, \
           'BLLAC'    : 3.19, \
           'MWC349A'  : 2.04, \
           'SGRA'     : 3.18 }
  return flux[src]
 
# --- convert day to day number --- #
def getDayNo( day ) :
  dayno = { '21mar' : '079', \
            '22mar' : '080', \
            '23mar' : '081', \
            '25mar' : '083', \
            '26mar' : '084', \
            '27mar' : '085', \
            '29mar' : '088', \
            '31mar' : '090', \
            '01apr' : '091', \
            '02apr' : '092', \
            '04apr' : '094' }
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
  
# --- generate calibration tables for all 5 days --- #
def doall2013() :

#  cpList = [2,3,4,5,6,8,9,14]
#  oneday( '21mar', cpList )
#  oneday( '22mar', cpList )

  cpList = [2,4,5,6,8,9,13,14]
#  oneday( '23mar', cpList )
#  oneday( '25mar', cpList )
  oneday( '26mar', cpList )
#  oneday( '27mar', cpList )

