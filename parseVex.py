# parseVex.py converts Haystack vex file to CARMA sched file

import time
import numpy

# search for these key words in vex file:
#  "scan" - contains scan number
#  "mode=1mmRDBE;" - contains source name, start time
#  "station=K1:" - contains duration
# don't process scan and mode lines until we see "station=K1", in case
#  CARMA is not involved in this scan


def UTtoLST( t ) :
  # t = time.mktime( [2013,03,13,01,07,30,0,0,0] ) 
    # allows for testing
  test = False 
  offset = 7958.
  t1 = t/.99726958
    # according to Wikipedia, each sidereal day is .99726958 mean solar days
    # sidereal time runs faster; find number of sidereal secs since Jan 1, 1970
  t1 = t1 - offset
    # offset was chosen to give correct answer on Mar 12, 2013
  if test :
    return time.strftime( "%H:%M:%S", time.localtime( t1 ) )
  else :
    return time.strftime( "%H:%M", time.localtime( t1 ) )


# makeSched converts .vex file to .Csched file
# reads entire vex file into lists, computes gaps and combines scans if needed,,
#   inserts pointing checks, etc, then writes scheduling file
# infile = name of vex file
# outfile = sched file
# intsecs = integration time per scan - putting it into the sched file allows us to
#      change this on the fly without restarting the script; default is 10 sec 

def makeSched( infile, outfile, intsecs=10 ) :
  fin = open(infile, 'r')
  fout = open(outfile, 'w')
  fout.write("#  %s created from file %s using parseVex.makeSched2\n#\n" % (outfile,infile) )
  fout.write("#  codes: vlb = VLBI scan\n")
  fout.write("#         reg = regular (calibration) scan\n" )
  fout.write("#         opt = optical pointing\n" )
  fout.write("#         pnt = radio pointing\n" )
  fout.write("#         xyp = xyphase calibraton\n" )
  fout.write("#         --- = nothing; go to next scan\n" )
  fout.write("#\n" )
  fout.write("#   source  UTstart    UTstop   int code   dur         scan         LST\n")

  scanNo = []
  start = []
  duration = []  
  source = []
  
  # --- loop through vex file, fill out lists with  scanNo, source, starttime --- #
  for line in fin:
    a = line.strip().split()
    if len(a) > 1 :
      
    # first line we stumble on will contain scan number
      if a[0] == "scan" :

      # before writing new scan, make sure that duration was filled in for previous scan;
      # if not, this indicates that CARMA was not used for previous scan - delete previous
      #   scan number, source, and start time from the lists
        if ( len(duration) < len(scanNo) ) :
          delScan = scanNo.pop()
          print "...deleting scan %s - CARMA not participating" % delScan
          start.pop()
          source.pop()

        scanNo.append( a[1][2:6] )  
        print "scan = %s" % scanNo[ len(scanNo)-1 ]

    # next line will contain source and starttime
    # save the start time as seconds since the epoch 
      if a[1].startswith("mode=") :
        start.append( time.mktime(time.strptime( a[0], "start=%Yy%jd%Hh%Mm%Ss;" ) ) )
        source.append( a[2][ 7 : len(a[2])-1 ] )
        #print "start time = %s" % ( time.strftime( "%H:%M:%S", time.localtime( start[ len(start)-1 ] )))
        #print "source = %s" % source[ len(source)-1 ]

    # last line will contain duration
      if a[0] == "station=Cm:" :
        duration.append( float(a[3]) ) 

  fin.close()

  # --- loop through lists 2nd time, compute time gap BEFORE each scan --- $
  print "\n=== BEGIN LOOP 2 ==="
  gap = numpy.zeros( [len(duration)], dtype=float )
  gap[0] = 3600.   # arbitrarily say that we start 60 min before the first scan
  n = 1
  while n < len(scanNo) :
    gap[n] = start[n] - (start[n-1] + duration[n-1])
    
  # if gap < 60 sec and source is identical to previous source, combine the scans
    if (gap[n] < 60.) and (source[n] == source[n-1]) :
      print "combining scans %s and %s" % ( scanNo[n-1], scanNo[n] )
      duration[n-1] = duration[n-1] + gap[n] + duration[n]   
      scanNo[n-1] = scanNo[n-1] + "/" + scanNo[n]
      scanNo.pop(n)
      start.pop(n)
      source.pop(n) 
      duration.pop(n)
    else :
      n = n + 1

  # --- loop through lists 3nd time, generate sched file --- #
  # now we have complete list of scans with sources, start times, durations, gaps
  # next, compute stop times and write out the sched file, inserting pointing placeholders in gaps

  print "\n=== BEGIN LOOP 3 ==="
  tlast_rpt = -1000.
  tlast_opt = -1000.
  tlast_xyp = -1000.
  code = "---"
  for n in range( 0, len(scanNo) ) :
    if (gap[n] > 10.*60.) and (start[n] > tlast_rpt + 2.*60.*60. ) : 
      code = "rpt"
      tlast_rpt = start[n]
    elif (gap[n] > 180.) and (start[n] > tlast_opt + 30.*60.)  :
      code = "opt"
      tlast_opt = start[n]
    elif (gap[n] > 180.) and (start[n] > tlast_xyp + 30.*60.) :
      code = "xyp"
      tlast_xyp = start[n]
    else :
      code = "---"
    stopstring = time.strftime( "%H:%M:%S", time.localtime( start[n] - 60.) )  
    fout.write("%10s     -      %s  %02ds  %3s  %5.2f\n" % \
	    (source[n], stopstring, intsecs, code, gap[n]/60. ))
    startstring = time.strftime( "%H:%M:%S", time.localtime( start[n] ) )
    stopstring = time.strftime( "%H:%M:%S", time.localtime( start[n] + duration[n] ) )
    code = "vlb"
    fout.write("%10s  %s  %s  %02ds  %3s  %5.2f   scan %9s  %s LST\n" % \
      ( source[n], startstring, stopstring, intsecs, code, duration[n]/60., \
      scanNo[n].ljust(9), UTtoLST( start[n] ) ) )
  fout.close()


# convert timeString to decimal hrs, without taking into account the day
def dechrs( timeString ) :
  [h1,m1,s1] = timeString.split(":")
  return int(h1) + float(m1)/60. + float(s1)/3600.


# makePlot shows schedule graphically
def makePlot( infile, outfile ) :
  color = { '0234+285' : 3, \
            '3C84'     : 4, \
            '3C111'    : 5, \
            '0721+713' : 6, \
            '0854+201' : 7, \
            'OJ287'    : 7, \
            '3C273'    : 8, \
            'M87'      : 2, \
            '3C279'    : 9, \
            '1633+382' : 3, \
            '3C345'    : 4, \
            'NRAO530'  : 5, \
            'SGRA'     : 2, \
            '1749+096' : 6, \
            '1921-293' : 7, \
            '2013+370' : 8, \
            'BLLAC'    : 9    }
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      if (len(a) > 5) and (a[4] == "vlb") :
        fout.write( "color %d\n" % (color[a[0]]) )
        fout.write( "fill 1\n" )
        fout.write( "rect %.3f %.3f -0.2 0.2\n" % (dechrs(a[1]),dechrs(a[2])) )   # start LST
        fout.write( "fill 2\n" )
        fout.write( "color 1\n" )
        fout.write( "rect %.3f %.3f -0.2 0.2\n" % (dechrs(a[1]),dechrs(a[2])) )   # start LST
  fin.close()
  fout.close()

  fout = open( "key.wip", "w" )
  xsize= 3
  xgap = 1
  ysize = .3
  ygap = .1

  srclist = color.keys() 
  x1 = 1
  y1 = 1
  for src in srclist :
    fout.write( "color %d\n" % (color[src]) )
    fout.write( "fill 1\n" )
    fout.write( "rect %.3f %.3f %.3f %.3f\n" % (x1,x1+xsize,y1,y1+ysize) )
    fout.write( "fill 2\n" )
    fout.write( "color 1\n" )
    fout.write( "rect %.3f %.3f %.3f %.3f\n" % (x1,x1+xsize,y1,y1+ysize) )
    fout.write( "move %.3f %.3f\n" % ( x1+xsize/2., y1+ysize/4. ) )
    fout.write( "putlabel 0.5 %s\n" % src)
    x1 = x1 + xsize + xgap
    if x1 > (23. - xsize) :
      x1 = 1
      y1 = y1 + ysize + ygap
  fout.close()

def done() :
  print "nothing"

def doit() :
  makeSched( "tr001.vex", "tr001.Csched", intsecs=10 ) 
  makeSched( "xx002.vex", "xx002.Csched", intsecs=10 ) 
  makeSched( "xx003.vex", "xx003.Csched", intsecs=10 ) 
  makeSched( "xx004.vex", "xx004.Csched", intsecs=10 ) 
  makeSched( "xx005.vex", "xx005.Csched", intsecs=10 ) 
  makePlot( "tr001.Csched", "tr001.wip" )
  makePlot( "xx002.Csched", "xx002.wip" )
  makePlot( "xx003.Csched", "xx003.wip" )
  makePlot( "xx004.Csched", "xx004.wip" )
  makePlot( "xx005.Csched", "xx005.wip" )


