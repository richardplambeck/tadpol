# flagGains.py
# produces file of uvflag commands to flag data where amplitude gains exceed some threshold

import numpy
import math
import subprocess
import shlex
import vCal
import time
import datetime
import matplotlib.pyplot as pyplot
import matplotlib.dates as mpdates

def getGains( visFile ) :
  '''produce time and gain arrays'''
  print "dump gain ratio from file %s" % visFile
  g = []         # accumulate gains in 1-D list of floats
  tt = []
  dayno = []
  type = None
  p= subprocess.Popen( ( shlex.split('gpplt vis=%s log=tmpGain' % visFile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  result = p.communicate()[0]
  print result
  fin = open( 'tmpGain', 'r' )
  for line in fin :
    if line.startswith( '# Listing of the Amp of the gains for I' ) :
      type = "I"
    if line.startswith( '# Listing of the Amp of the gains for X,Y' ) :
      type = "XY"
    if (len(line) > 0) and (not line.startswith("#")) :
      a = line.split()
      if len(a) == 8 :  # begins new record
        dayno.append( int(a[0]) )
        tt.append( a[1] )
        del a[0:2]      # discard daynumber and UT
      for i in range(0,len(a)) :
        g.append( float( a[i] ) )
  fin.close()
  #if type == "XY" :
  #  gg = numpy.reshape( numpy.array(g), [len(tt),-1, 2] )     # gg[time, [ant, [GR,GL]]]
  gg = numpy.reshape( numpy.array(g), [len(tt),-1] )     # gg[times][ants]
  if type == "I" :
    ggg = numpy.delete( gg, numpy.s_[15:], axis=1 )        # trim off gains for ants > 15
  else :
    ggg = numpy.delete( gg, numpy.s_[30:], axis=1 )        # trim off gains for ants > 15
  return [tt, ggg]

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
  [t, g ] = getGains( infile ) 
  type = "I"
  if len(numpy.shape(g)) > 2 :
    type = "XY"
  fout = open( "flag.csh", "w" )
  fout.write ("#!/bin/csh\n\n")
  fout.write ("set FILE = %s\n\n" % infile)
  for tt,gg in zip( t, g) :
    gstr = ""
    for n,ggg in enumerate(gg) :
      if (len(gg) > 15) and (n % 2 == 0) :
        gstr = gstr + " "
      gstr = gstr + "%5.1f" % ggg
    bad = " "
    if gg.max() > cutoff :
      bad = "*"
    
      t1 = t2 = tt
      fout.write("uvflag flagval=flag vis=$FILE select='%s'\n" % (timeString(t1) ) )
      fout.write("   # %s\n" % gstr )
    print "%s %s %s" % (tt,bad,gstr)
  fout.close()

antList0 = [1, 2,3,4,5,6,7,8,9,11,12,13,15]
# --- read varplt log containing tsys or elev, 6 lines per vis record, 4 values per line --- #
def flagTsys( infile, cutoff=2500., antList=antList0 ) :
  varlist = []
  t = []
  fin = open( infile, "r" )
  n1 = 0
  nline = 0
  for line in fin:
    if ( not line.startswith("#") ) :
      a = line.split()
      if ( len(a) == 6) :
        nline = nline + 1
        t.append(datetime.datetime.strptime( a[1], "%H:%M:%S") )
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
  vl = numpy.reshape( numpy.array( varlist, dtype=float ), (-1,15) )
  va = numpy.hsplit( vl, 15 )
  #print vl[13][0:]
  print va[13]
  plot( t, va, cutoff, antList )

def plot( t, s, cutoff, antList) :
  color = [ "red", "blue", "green", "chartreuse", "orangered", \
            "aqua", "fuchsia", "gray", "lime", "maroon", "navy", \
            "olive", "orange", "silver", "teal", "black" ]
  pyplot.clf()
  pyplot.ion()
  fig = pyplot.subplot(1,1,1)
  tt = mpdates.date2num(t)
  pyplot.ylim(0,10000.)
  for ant,col in zip(antList,color) :
    fig.plot_date( t, s[ant-1], ".", color=col, label="C%d" % ant )
    for ttt,tsys in zip( t, s[ant-1] ) :
      if tsys > cutoff :
        print "   C%d   %s   %.0f" % (ant, ttt.strftime( "%H:%M:%S" ), tsys )
      
  fig.legend( loc=0, prop={'size':10} )

  #fig.xaxis.set_major_locator(mpdates.MinuteLocator(interval=10))
  #fig.xaxis.set_major_formatter(mpdates.DateFormatter( "%H:%M" ) )
  fig.grid()
  pyplot.ylabel("Tsys (K)")
  pyplot.xlabel("UT")
  #fig.xaxis.set_major_locator(mpdates.AutoDateLocator())
  #fig.xaxis.set_major_formatter(mpdates.AutoDateFormatter())
  pyplot.show()

#flagTsys( "win1LL.log" )
#flagTsys( "win1RR.log" )
