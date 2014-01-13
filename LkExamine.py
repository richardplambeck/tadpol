# LkExamine.py
# goes through a leakage file
# prints sum of DR and DL for each freq interval
# prints avg and median magnitudes for each freq interval

import numpy
#import sys

def mags( oneLk ) :
  fin = open( oneLk )
  DRsum = 0.+0.j
  DLsum = 0.+0.j
  DRmaglist = []
  DLmaglist = []
  for line in fin :
    if not line.startswith("#") :
      a = line.split()
      if len(a) > 7 :
        DR = float(a[3]) + float(a[4])*1j
        DL = float(a[5]) + float(a[6])*1j
        DRmag = abs(DR)
        DLmag = abs(DL)
        if (DRmag != 0.) or (DLmag != 0.) :
          DRmaglist.append(DRmag)
          DLmaglist.append(DLmag)
          DRsum = DRsum + DR
          DLsum = DLsum + DL
        print "  %s  (%+5.3f%+5.3fj)  (%+5.3f%+5.3fj)  %5.3f  %5.3f" % ( a[0], DR.real, DR.imag, DL.real, DL.imag, DRmag, DLmag )
        if a[0] == "C15" :
          #print " sums       ", DRsum, DLsum
          print " sums  (%+5.3f%+5.3fj)  (%+5.3f%+5.3fj)" % ( DRsum.real, DRsum.imag, DLsum.real, DLsum.imag )
          print " mean mag     %5.3f  %5.3f" % ( numpy.mean(DRmaglist),numpy.mean(DLmaglist) )
          print " median mag   %5.3f  %5.3f" % ( numpy.median(DRmaglist),numpy.median(DLmaglist) )
          DRsum = 0.+0.j
          DLsum = 0.+0.j
          DRmaglist = []
          DLmaglist = []
  fin.close()

# copy leakages from inLk to outLk, referencing all leakages to DR(zant),DL(zant)
def zero( inLk, outLk, zant ) :
  fin = open( inLk, "r" )
  fout = open( outLk, "w" )
  DR = numpy.zeros( 16, dtype=complex )
  DL = numpy.zeros( 16, dtype=complex )
  for line in fin :
    if line.startswith("#") :
	  fout.write( "%s" % line )
    else :
      a = line.split()
      if (len(a) > 9) and (line.startswith("C")) :        # new style Lk table, includes all antennas
        ant = int( a[0].strip("C") )
        f1 = min(float(a[1]),float(a[2]))
        f2 = max(float(a[1]),float(a[2]))
        DR[ant] = float(a[3]) + 1j * float(a[4])
        DL[ant] = float(a[5]) + 1j * float(a[6])
        lineString = a[9]
        if ant == 15 :
          for n in range(1,16) :
            if (abs(DR[n]) > 0) or (abs(DL[n]) > 0) :   # leave zeros untouched!
              DR[n] = DR[n] - DR[zant] + .001           # add .001 to distinguish from zero!
              DL[n] = DL[n] - DL[zant] + .001
            fout.write("C%02d %8.3f %8.3f %8.3f %6.3f %8.3f %6.3f  %s  %s    %s\n" % \
              ( n, f1, f2, DR[n].real, DR[n].imag, DL[n].real, DL[n].imag, a[7], a[8], a[9]) )
  fin.close()
  fout.close()
 
     

