# scuvec.py
# reads SCUBA data file (table6.txt), writes olay file for cgdisp plots

import math

srcname = "Serpens"
radius = 10.

def scuba() :
  scale = 2.   # 2% pol = 1 arcsec on plot
  fin = open("table6.txt", "r")
  fout = open( "SCUBAOlay", "w" )
  nline = 0
  for line in fin:
    nline = nline + 1
    # if line.startswith("Serpens") :
    if True :
      print line
      RA = line[33:45]
      DEC = line[45:57]
      polpct = float(line[77:82])
      PA = float(line[86:92])
      PAerr = float(line[92:])
      if PAerr < 15. :
        PA = PA + 90.   # rotate to show B-field direction
        print "RA =",RA 
        print "DEC =",DEC
        print "polpct = %.1f" % polpct
        print "PA = %.1f\n" % PA
        #fout.write("ocircle hms dms %d no %s %s %.1f\n" % (nline,RA,DEC,radius) )
        [ra1,dec1,ra2,dec2] = vec(RA,DEC,polpct,PA,scale)
        fout.write("line hms dms %d no %s %s %s %s\n" % (nline,ra1,dec1,ra2,dec2) )
  fin.close()
  fout.close()

# vec computes ra and dec of line segment ends for pol vector
# scale is arcsec per 1% pol

def vec(RA,DEC,polpct,PA,scale) :
  a = RA.split()
  rhr = float(a[0])
  rmin = float(a[1])
  rsec = float(a[2])  
  print rhr, rmin, rsec
  a = DEC.split()
  ddeg = float(a[0])
  dmin = float(a[1])
  dsec = float(a[2]) 
  if (ddeg < 0.) :
    dmin = -1.*dmin
    dsec = -1.*dsec
  #print ddeg, dmin, dsec

  # RA and DEC in radians
  RArad = math.pi * (rhr + rmin/60. + rsec/3600.)/12.
  DECrad = math.pi * (ddeg + dmin/60. + dsec/3600.)/180.
  
  # length of line segment in arcsec
  length = scale * polpct * math.pi/(180.*3600.)
  dra = math.sin(math.pi*PA/180.)/math.cos(DECrad) * length/2.
  ddec = math.cos(math.pi*PA/180.) * length/2. 
  #print dra, ddec
  return [ hms(RArad-dra), dms(DECrad-ddec), hms(RArad+dra), dms(DECrad+ddec) ] 

# convert radians into hh:mm:ss.ss
def hms(rad) :
  hh = int(12.*rad/math.pi)
  rad = rad - math.pi*hh/12.
  min = int(60.*12.*rad/math.pi)
  rad = rad - math.pi*min/(12.*60.)
  sec = 12.*3600.*rad/math.pi
  return " %2d %2d %.3f" % (hh,min,sec)

def dms(rad) :
  deg = int(180.*rad/math.pi)
  rad = rad - math.pi*deg/180.
  min = int(60.*180.*rad/math.pi)
  rad = rad - math.pi*min/(60.*180.)
  sec = 3600.*180.*rad/math.pi
  if min < 0. :
    min = -1.*min
  if sec < 0. :
    sec = -1.*sec
  return " %3d %2d %.2f" % (deg,min,sec)

def hertz() :
  scale = 4.   # 1% pol corresponds to 4 arcsec on map

  # create dictionary of ra,dec vs src name
  raSrc = dict()
  decSrc = dict()
  f1 = open("/fringe2/plambeck/pol/Hertz/table1.dat")
  for line in f1 :
    src = line[0:15]
    rhr = float(line[54:56])
    rmin = float(line[57:59])
    rsec = float(line[60:64])
    ddeg = float(line[65:68])
    dmin = float(line[69:71])
    dsec = float(line[72:74])
    if ddeg < 0 :
      dmin = -1. * dmin
      dsec = -1. * dsec
    print "%s  %3.0f %2.0f %4.2f  %3.0f %3.0f %4.1f" % (src,rhr,rmin,rsec,ddeg,dmin,dsec) 
    raSrc[ src ] = math.pi * (rhr + rmin/60. + rsec/3600.)/12.
    decSrc[ src ] = math.pi * (ddeg + dmin/60. + dsec/3600.)/180.
  f1.close()
  
  # now read through table2, generate olay lines
  nline = 0
  f1 = open("/fringe2/plambeck/pol/Hertz/table2.dat")
  fout = open("HertzOlay", "w")
  for line in f1:
    nline = nline + 1
    src = line[0:15]
    draArcsec = float(line[16:20])
    ddecArcsec = float(line[21:25])
    polpct = float(line[38:43])
    PA = float(line[50:55])
    PAerr = float(line[56:62])
    if PAerr < 15. :
      decRadians = decSrc[ src ] + math.pi * ddecArcsec/(180.*3600.) 
      raRadians = raSrc[ src ] + math.pi * draArcsec / (math.cos(decRadians) * 180.*3600.)
      RA = hms(raRadians)
      DEC = dms(decRadians)
      PA = PA + 90.   # rotate to show B-field direction
      print "RA =",RA 
      print "DEC =",DEC
      print "polpct = %.1f" % polpct
      print "PA = %.1f\n" % PA
      #fout.write("ocircle hms dms %d no %s %s %.1f\n" % (nline,RA,DEC,radius) )
      [ra1,dec1,ra2,dec2] = vec(RA,DEC,polpct,PA,scale)
      fout.write("line hms dms %d no %s %s %s %s\n" % (nline,ra1,dec1,ra2,dec2) )
  f1.close()
    
    
 
	
