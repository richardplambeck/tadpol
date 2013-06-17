import math
import time
import numpy
import subprocess
import shlex

result = 0

def dechrs( hhmmss ) :
  hr,min,sec = hhmmss.split(":")
  return float(hr) + float(min)/60. + float(sec)/3600.

def hhmmss( dechrs ) :
  hrs = int(dechrs)
  decmins = 60.*(dechrs - hrs)
  mins = int(decmins)
  decsecs = 60.*(decmins - mins)
  secs = int(decsecs)
  return "%02d:%02d:%02d" % (hrs,mins,secs)

def old() :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  actual = 0
  target = 0
  for line in fin :
    a = line.split()
    if (len(a) >= 10) and (a[8] == "vlbi") :
      targetStart = a[3]
      #print targetStart
      target = dechrs( targetStart )
      targetSource = a[1]
      actual = 0
    if (len(a) >= 7) and (a[3] == "vlbi") :
      actualStart =  a[0][8:] 
      #print actualStart
      actual = dechrs( actualStart )
      if (actual > target ) :
        print "target ",targetStart,"   actual ",actualStart, "   ",targetSource
        fout.write("target %s   actual %s   source %s\n" % (targetStart,actualStart,targetSource) )
      else :
        print "OK"
        fout.write("OK\n")
  fin.close()
  fout.close()

def schedScans( schedfile ) :
  schedList = []
  fin = open( schedfile, "r" )
  for line in fin :
     if len(line) > 1 :
       a = line.split()
       if len(a) > 8 and a[6] == "scan" :
         print a[7]
         scan = { "number" : a[7], "targetStart" : a[1], "source" : a[0] }
         schedList.append( scan )
  fin.close
  return schedList

# reads through the entire vlbiLog to search for scan with same wanted
#    start time, with format "yymmmdd:hh:mm:ss")
# returns actual start timein same format
# if scan is not found, return "missed"
def findScan( wantedstart, logfile="vlbiLog" ) :
  fin = open( logfile, "r" )
  for line in fin :
    if len(line) > 1 :
      a = line.split()
      if (len(a) > 9) and (a[8] == "vlbi") : 
        schedstart = a[0][0:8] + a[3]
        #print schedstart
      elif (len(a) > 4) and (a[3] == "vlbi") :
        actualstart = a[0]
        if schedstart == wantedstart :
          return actualstart   
  return "missed"   # if we drop through the file without finding this scan

# findLate finds late or missed scans in a schedule file
def findLate( date, schedfile, outfile ) :
  fout = open( outfile, "w" ) 
  fout.write("# look for missing or late scans in %s\n" % schedfile )
  schedList = schedScans( schedfile )
  for scan in schedList :
    schedstart = date + scan["targetStart"]
    actualstart = findScan( schedstart )
    print "actualstart = ", actualstart
    print "%10s %8s  %s  %s" % ( scan["number"], scan["source"], scan["targetStart"], actualstart)      
    if actualstart == "missed" :
       fout.write("%10s  %8s  %s  %s\n" % ( scan["number"], scan["source"], scan["targetStart"], actualstart) )      
    else :
       target = dechrs( scan["targetStart"] )
       actual = dechrs( actualstart[8:] )
       fout.write("%10s  %8s  %s  %s   %6.1f" % ( scan["number"], scan["source"],  \
           scan["targetStart"], actualstart[8:], 3600.*(actual-target) ) )
       if target > actual :
         fout.write("  OK")
       else :
         fout.write("  secs LATE")
       [ p90, p180, frate ] = phaseSwitchState( actualstart )
       if p90 == -1 :
         fout.write( " NO DATABASE ENTRY  ")
       else :
         if (p90 != 1) : fout.write( "  P90 on  ")
         if (p180 != 1) : fout.write( "  P180 on  ")
         if (frate != 114.) : fout.write( "  no 114 Hz  ")
       fout.write("\n")
  fout.close() 

# queries the MPstore halfsecond database to learn if 90 and 180 degree phase
#  switching wree turned off, and a 114 Hz signal added, to C1 during vlbi scan
# actualstart has format "13mar23:12:01:50"

def phaseSwitchState( actualstart ) :
  yr = "20%2s" % actualstart[0:2]
  mon = actualstart[2:5]
  day = actualstart[5:7]   
  tstart = dechrs( actualstart[8:] ) 
  t1 = "%s-%s-%s %s" % (yr, mon, day, hhmmss(tstart + 10./3600.) )   # 10 sec after start
  t2 = "%s-%s-%s %s" % (yr, mon, day, hhmmss(tstart + 12./3600.) )  # 12 sec after start
  print " checking MPstore  ", t1, t2
  pastHeader = False
  p90 = -1
  p180 = -1
  frate = -1.
  
  p= subprocess.Popen( shlex.split('/home/eml/bin/readArcNew.csh \
"array.frame.utc carmastring;array.frame.validity;\
Loberotator.Channel1.phaseSwitch90;\
Loberotator.Channel1.phaseSwitch180;\
Loberotator.Channel1.offsetPhaseRate" \
"%s" "%s" 0' % (t1, t2) ), stdout=subprocess.PIPE, stderr=subprocess.PIPE )
  result = p.communicate()[0]
  try :
    lines = result.split("\n")
    for line in lines :
      if len(line) > 0 :
        a = line.split()
        if (len(a) > 5) and pastHeader :
          p90 = int(a[2])
          p180 = int(a[3])
          frate = float(a[4])
          return [p90, p180, frate]
        elif (len(a) > 1) and (a[0] == "Date/Time") :
          pastHeader = True 
    print "error querying database"
    return [-1, -1, -1]
  except :
    print "error querying database"
    return[ -1, -1 ,-1 ]   
    

def doit() :
  #findLate( "13mar21:", "21mar/Sched.21mar", "21mar/lateList.txt" )
  #findLate( "13mar22:", "22mar/Sched.22mar", "22mar/lateList.txt" )
  #findLate( "13mar23:", "23mar/Sched.23mar", "23mar/lateList.txt" )
  #findLate( "13mar25:", "25mar/Sched.25mar", "25mar/lateList.txt" )
  #findLate( "13mar26:", "26mar/Sched.26mar", "26mar/lateList.txt" )
  findLate( "13mar27:", "27mar/Sched.27mar", "27mar/lateList.txt" )

def didit() :
  list = phaseSwitchState( "13mar26:15:46:39" )
  print list
