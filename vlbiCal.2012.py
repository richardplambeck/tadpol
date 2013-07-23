# vlbiCal.py
#
# computes SEFD for phased and unphased antennas during vlbi run
# 2012 version - deleted emGains...

import math
import time
import numpy
import subprocess
import shlex

# --- ants phased in 2012
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
        #print a[1], t[nline]
        if ( (nline > len(t)) or (abs(dechrs(a[1]) - t[nline]) > .0003)) :
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

# --- wrapper for avgSEFD that uses default gains if no pre- or post-scan data --- #
def scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsys, antList ) :
  s = avgSEFD( t1, t2, src, t, source, gain, tsys, antList ) 
  if (s > 0.) :
    fout.write("%8.0f " % (100.*round(s/100.)) )
  else :
    s = avgSEFD( t1, t2, src, t, source, defgain, tsys, antList )  
    fout.write("%8.0f*" % (100.*round(s/100.)) )
  return s

# --- compute SEFD averaged over a scan --- #
# ... t1 and t2 are decimal hrs since UT0, src is source name
# ... t is (npts x 1) time array, decimal hrs since UT0
# ... source is (npts x 1) list of source names
# ... gain is (npts x 15) complex array
# ... tsys is (npts x 15) real array
def avgSEFD( t1, t2, src, t, source, gain, tsys, antList ) : 
  print "t1,t2,src,antList = %.4f,%.4f,%s,%s" % (t1,t2,src,antList)
  navg = 0
  avg = 0.
  for n in range( 0, len(t) ) :
    if ( (t[n] >= t1) and (t[n] <= t2) and (src == source[n]) ) :
      s = SEFD( gain[n], tsys[n], antList )
        # ... note that gain[n] and tsys[n] are (1 x 15) arrays
      print "%8.4f %6.3f %6.0f %9.0f" % (t[n], abs(gain[n][0]), tsys[n][0], s)
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
 

# --- generate calibration tables for all 5 days --- #
def doall() :
  
  cpList = [2,3,5,9,11,12,13]
  oneday( '04apr', cpList )
  emoneday( '04apr', cpList )

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


# --- generate calibrations for one day (e.g., day='day075')--- #
# ... cpList is list of antennas that are phased (used if summ file says 'CP')
# ... tsysL0,tsysR0 are the 'lo-band' systemps (from downconverters 5 and 6)
# ... tsysL1,tsysR1 are the 'hi-band' systemps (from downconverters 7 and 8)

def oneday( day, cpList=antList ) :

  # ... retrieve time and gain arrays
  [time,gain] = getGains( day+"/wide.cal" )

  # ... create decimal time array 
  t = numpy.empty( len(time), dtype=float )
  for n in range (0, len(time) ) :
    t[n] = dechrs( time[n] )

  # ... retrieve tsys, rmspath, tau230, elev arrays
  tsysL0 = getVar15( day+"/tsys13L.log", t )
  tsysR0 = getVar15( day+"/tsys13R.log", t )
  tsysL1 = getVar15( day+"/tsys15L.log", t )
  tsysR1 = getVar15( day+"/tsys15R.log", t )
  rmspath = getVar( day+"/rmspath", t )
  tau230 = getVar( day+"/tau230", t )
  source = getSource ( day+"/uvindex.log", t )
  elevarray = getVar15( day+"/elev", t )

  # ... use elev of C1 as the elevation
  elev = numpy.empty( (len(time)), dtype=float )
  for n in range (0, len(time) ) :
    elev[n] = elevarray[n][0]

  # ... all array lengths should match, otherwise this routine will crash!
  print len(time), len(gain), len(tsysL0), len(tsysR0), len(tsysL1), \
    len(tsysR1), len(rmspath), len(tau230), len(source), len(elev)

  # ... create absgain and default gain files for future use
  absgain = numpy.empty( (len(time),15), dtype=complex )
  defgain = numpy.empty( (len(time),15), dtype=complex )
  jyperk = numpy.array( [1.07,1.12,1.02,1.09,1.09,1.30,1.03,1.06,1.03,1.,1.34,1.09,1.04,1.03,1.05] )
  for n in range (0, len(time)) : 
    for i in range (0,14) :
      absgain[n][i] = numpy.abs( gain[n][i] ) + 1j * 0.
      defgain[n][i] = jyperk[i] * (1. + 0j)
  print "\ndefgains:\n",defgain

  # ... make table of gain vs elev for plotting
  gainVsElev( elevarray, gain, day+".gainVsElev" )

  # ... process the 'summ' file which lists the schedule
  fin = open( day+"/summ", "r" )
  fout = open( day+".results", "w" )
  fout.write("#   scan     No       src   UTstart    UTstop  el  tau  path   CL-lo    FL-lo    FL-hi    FR-lo    FR-hi   pheff\n")
  fout.write("#                                             deg        um   SEFD-Jy  SEFD-Jy  SEFD-Jy  SEFD-Jy  SEFD-Jy\n")
  fout.write("#\n")
  for line in fin :

    # ... print scanname, src, UTstart, UTstop (from summ file)
    if not line.startswith("#") :
      a = line.split()
      src = a[0]
      utStart = a[1]
      utStop = a[2]
      
      scanname = a[9]
      scanNo = a[10][0:4]
      print scanname, src
      fout.write(" %8s  %4s %9s  %8s  %8s " % (scanname,scanNo,src,utStart,utStop) ) 

      # ... create decimal start and stop times, print tau230, rmspath
      t1 = dechrs( utStart )
      t2 = dechrs( utStop )
      elevavg = avgVar( t1, t2, t, elev )
      tau230avg = avgVar( t1, t2, t, tau230 )
      rmspathavg = avgVar( t1, t2, t, rmspath )
      fout.write(" %2d %5.2f %4.0f" % (elevavg,tau230avg,rmspathavg ) ) 


      # ... print SEFD for C1 lo band only (hi not recorded)
      print "\nC1"
      s10 = scanSEFD( fout, t1-.1, t2+.1, src, t, source, defgain, defgain, tsysL0, [1] ) 

      # ... print SEFDs for phased array L and phased array R, lo and hi bands
      sL0 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysL0, cpList ) 
      sL0_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysL0, cpList ) 

      sL1 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysL1, cpList ) 
      #sL1_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysL1, cpList ) 

      sR0 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysR0, cpList ) 
      #sR0_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysR0, cpList ) 

      sR1 = scanSEFD( fout, t1, t2, src, t, source, gain, defgain, tsysR1, cpList ) 
      #sR1_ideal = avgSEFD( t1, t2, src, t, source, absgain, tsysR1, cpList ) 

    # all efficiencies are the same, so only write out the L-lo one
      fout.write("   %4.2f "  % efficiency( sL0, sL0_ideal ) )

      #... print predicted lo- and hi-band CARMA-CARMA correlations
      #if (s10*sL0) > 0.1 :
      #  amplo = 1.e4 * sourceFlux(src) / math.sqrt(s10*sL0)
      #else :
      #  amplo = 0.
      #fout.write("  %6.3f" % amplo )

      fout.write("\n")

  # ... finished
  fin.close()
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

def doit( ) :
  basicEffic( antList, "day075/s.vis", "day075/s.phEffic" )
  cumEffic( 3, "day075/s.phEffic", "day075cum.wip" ) 
  basicEffic( antList, "day080/s.vis", "day080/s.phEffic" )
  cumEffic( 3, "day080/s.phEffic", "day080cum.wip" ) 
  basicEffic( antList, "day081/s.vis", "day081/s.phEffic" )
  cumEffic( 3, "day081/s.phEffic", "day081cum.wip" ) 

  
