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
        DRmaglist.append(DRmag)
        DLmaglist.append(DLmag)
        DRsum = DRsum + float(a[3]) + float(a[4])*1j
        DLsum = DLsum + float(a[5]) + float(a[6])*1j
        print DR,DL,DRmag,DLmag
        if a[0] == "C15" :
          print " sums       ", DRsum, DLsum
          print " mean mag   ", numpy.mean(DRmaglist),numpy.mean(DLmaglist)
          print " median mag ", numpy.median(DRmaglist),numpy.median(DLmaglist)
          DRsum = 0.+0.j
          DLsum = 0.+0.j
          DRmaglist = []
          DLmaglist = []
  fin.close()

