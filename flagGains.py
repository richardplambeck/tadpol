# flagGains.py
# produces file of uvflag commands to flag data where amplitude gains exceed some threshold

import numpy
import math
import vlbiCal

def flag( infile, cutoff ) :
  [time, gainComplex] = vlbiCal.getGains( infile ) 
  gainMag = numpy.abs( gainComplex )
  for t,g in zip( time, gainMag ) :
    if g.max() > cutoff :
      s = "[ "
      for n in range(0,15) :
        if g[n] > cutoff : s = s + "*"
        else : s = s + " "
      s = s + "]"
      print t, s 

# turns "hh:mm:ss" into "hh:mm:ss,HH:MM:SS" where t2 = t1 + 1 sec
def timeString( tstr ) :
  [strhr,strmin,strsec] = tstr.split(":")
  hr = int(strhr)
  min = int(strmin)
  sec = float(strsec) - 1.
  if sec < 0. :
    sec = sec + 60.
    min = min - 1
  if min < 0 :
    min = min + 60
    hr = hr - 1
  if (hr < 0) :
    print "### WARNING: WRAPPED THROUGH PREVIOUS DAY ###"
    return [ ]
  t1str = "%0.2d:%0.2d:%02d" % (hr,min,sec)
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
  t2str = "%0.2d:%0.2d:%02d" % (hr,min,sec)
  return "time(%s,%s)" % (t1str, t2str) 

def makeFile ( infile, cutoff ) :
  [time, gainComplex] = vlbiCal.getGains( infile ) 
  fout = open( "flag.csh", "w" )
  fout.write ("#!/bin/csh\n\n")
  fout.write ("set FILE = %s\n\n" % infile)
  gainMag = numpy.abs( gainComplex )
  for t,g in zip( time, gainMag ) :
    if g.max() > cutoff :
      t1 = t2 = t
      fout.write("uvflag flagval=flag vis=$FILE select='%s'\n" % (timeString(t1) ) )
      fout.write("   # %s\n" % (numpy.array_str( g, precision=2, max_line_width=200) ) )
  fout.close()
