# flagGains.py
# produces file of uvflag commands to flag data where amplitude gains exceed some threshold

import numpy
import math
import vlbiCal
import subprocess
import shlex

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
