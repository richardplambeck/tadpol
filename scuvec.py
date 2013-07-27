# scuvec.py
# reads SCUBA data file (table6.txt), writes olay file for cgdisp plots

import math

srcname = "Serpens"
radius = 10.

def scuba() :
  scale = 1.5   # 1% pol = 2 arcsec on plot
  fin = open("table6.txt", "r")
  fout = open( "SCUBAOlay", "w" )
  fout2 = open( "SCUBA.wip", "w" )
  nline = 0
  for line in fin:
    nline = nline + 1
    # if line.startswith("Serpens") :
    if True :
      RA = line[33:45]
      DEC = line[45:57]
      polpct = float(line[77:82])

      # quickie way of plotting all vectors the same length
      # polpct = 4
      # quickie way of plotting all vectors the same length

      PA = float(line[86:92])
      PAerr = float(line[92:])
      if PAerr < 15. :
        PA = PA + 90.   # rotate to show B-field direction
        print "RA,DEC,polpct,PA,PAerr = %s %s %.1f %.1f %.1f" % (RA,DEC,polpct,PA,PAerr)
        [ra1,dec1,ra2,dec2] = vec(RA,DEC,polpct,PA,scale)
        fout.write("line hms dms %d no %s %s %s %s\n" % (nline,ra1,dec1,ra2,dec2) )
        if (PAerr > 10.) :
          fout2.write("lwidth 1\n")
        fout2.write("radec %s %s\n" % (ra1,dec1))
        fout2.write("drawradec %s %s\n" % (ra2,dec2))
        if (PAerr > 10.) :
          fout2.write("lwidth 5\n")
  fin.close()
  fout.close()
  fout2.close()

# vec computes ra and dec of line segment ends for pol vector
# scale is arcsec per 1% pol

def vec(RA,DEC,polpct,PA,scale) :
  a = RA.split()
  rhr = float(a[0])
  rmin = float(a[1])
  rsec = float(a[2])  
  #print rhr, rmin, rsec
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

# output table in form that Chat wants: ra_degr, dec_degr, pct, pct_err, PA_deg, PA_err
def hertz2() :
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
  fout = open("Hertz.dat", "w")
  for line in f1:
    src = line[0:15]
    draArcsec = float(line[16:20])
    ddecArcsec = float(line[21:25])
    polpct = float(line[38:43])
    polpcterr = float(line[43:49])
    PA = float(line[50:55])
    PAerr = float(line[56:62])
    if PAerr < 15. :
      decRadians = decSrc[ src ] + math.pi * ddecArcsec/(180.*3600.) 
      raRadians = raSrc[ src ] + math.pi * draArcsec / (math.cos(decRadians) * 180.*3600.)
      RA = hms(raRadians)
      RAdeg = 180./math.pi * raRadians
      DEC = dms(decRadians)
      DECdeg = 180./math.pi * decRadians

      PA = PA + 90.   # rotate to show B-field direction
      print "RA =",RA 
      print "DEC =",DEC
      print "polpct = %.1f" % polpct
      print "PA = %.1f\n" % PA
      fout.write("%10.5f %10.5f  %5.2f %5.2f %8.1f %5.1f     # %15s %4.0f %4.0f %s %s\n" % 
        (RAdeg, DECdeg, polpct, polpcterr, PA, PAerr, src, draArcsec, ddecArcsec, RA, DEC) )
  f1.close()
  fout.close()


# special routine to process SHARP data from Davidsen et al table 1
def L1527() :
  scale = 3.   # 1% pol corresponds to 4 arcsec on map
  rhr = 4.
  rmin = 39.
  rsec = 53.9
  ddeg = 26.
  dmin = 3.
  dsec = 11.
  raSrc = math.pi * (rhr + rmin/60. + rsec/3600.)/12.
  decSrc = math.pi * (ddeg + dmin/60. + dsec/3600.)/180.
  
  nline = 0
  f1 = open("/fringe2/plambeck/pol/Hertz/L1527SHARP.dat")
  fout = open("SHARPOlay", "w")
  for line in f1:
    a = line.split()
    print a
    draArcsec = float(a[0])
    ddecArcsec = float(a[1])
    polpct = float(a[2])
    PA = float(a[4])
    PAerr = float(a[5])
    if PAerr < 15. :
      decRadians = decSrc + math.pi * ddecArcsec/(180.*3600.) 
      raRadians = raSrc + math.pi * draArcsec / (math.cos(decRadians) * 180.*3600.)
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
  fout.close()
    
# generate arcmin grid, staggering az or el

def gridMake( delArcmin, staggerEl=False ) :
  fout = open("grid", "w")
  deltaAz = delArcmin
  deltaEl = 0.866 * delArcmin
  if staggerEl :
    deltaAz = 0.866 * delArcmin
    deltaEl = delArcmin 
  for m in range( -10, 11, 1 ) :
    for n in range( -10, 11, 1 ) :
      off = 0.
      daz = n * deltaAz
      dele = m * deltaEl
      if staggerEl and (n % 2 == 1) :
        dele = dele + 0.5 * deltaEl
      elif (m % 2 == 1) :
        print "staggering this row"
        daz = daz + 0.5 * deltaAz
      print daz, dele
      fout.write("%9.3f %9.3f\n" % (daz,dele)  )
  fout.close()



def gridPlot( ) :
  ramax = math.pi * (18. + 30./60. + 0.7/3600.)/12.
  ramin = math.pi * (18. + 29./60. + 55.7/3600.)/12.
  decmin = math.pi * ( 1. + 12./60. + 50./3600.)/180.
  decmax = math.pi * ( 1. + 14./60. + 40./3600.)/180.
  print ramin,ramax,decmin,decmax
  ra0 = (ramax + ramin) / 2.
  dec0 = (decmin + decmax) / 2.
  print "center: ",hms(ra0),dms(dec0)
  fin = open("grid", "r")
  fout = open("gridOlay", "w")
  fout2 = open("gridList", "w")
  for line in fin :
    a = line.split()
    ra = ra0 + (math.pi * float(a[0]) /60 /12./15.) / math.cos(dec0)
    dec = dec0 + math.pi*float(a[1])/60./180.
    print a[0],a[1],hms(ra),dms(dec) 
    if (ra >= ramin) and (ra <= ramax) and (dec >= decmin) and (dec <= decmax) :
      print ">>",a[0],a[1],ra,dec
      fout.write("ocircle hms dms no no %s %s 15.\n" % (hms(ra),dms(dec))) 
      fout2.write("%s" % line)
  fin.close()
  fout.close()
  fout2.close()
 
# special routine to process MMS6 data from Matthews et al. 2005, Table1	
def mms6():
  ra0 = math.pi * (5. + 35./60. + 23.498/3600.)/12.
  dec0 = math.pi * (-5.  -1./60. - 32.221/3600.)/180.
  scale = 0.2
  fin = open( "/fringe2/plambeck/pol/MMS6/MMS6bima.dat", "r")
  fout = open( "MMS6.olay", "w" )
  nline = 0
  for line in fin:
    a = line.split()
    nline = nline + 1
    draArcsec = float(a[0])
    ddecArcsec = float(a[1]) 
    polpct = float(a[2])
    PA = float(a[4])
    decRadians = dec0 + math.pi * ddecArcsec/(180.*3600.) 
    raRadians = ra0 + math.pi * draArcsec / (math.cos(decRadians) * 180.*3600.)
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
  fin.close()  
  fout.close()
 
