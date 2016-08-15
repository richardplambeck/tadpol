# ob.py
# these are routines for dealing with the Optics Bench scans

import numpy
import math
import sys
import string
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages


def findminima( x, p, fGHz ):
    delta = 0.45*300./fGHz
      # anticipated spacing between minima
    xminlist = []
    pminlist = []
    while numpy.ma.count( x ) > 0 :
      nmin = numpy.argmin( p )
        # nmin is index of point with  minimum power
      xmin = x[nmin]
      xminlist.append( xmin )
      pminlist.append( p[nmin] )

    # now mask all points that are within +/-delta of the min
      x2 = numpy.ma.masked_inside( x, xmin-delta, xmin+delta )
      p2 = numpy.ma.masked_where( numpy.ma.getmask(x2), p )

    # copy the arrays to allow easy iteration
      x = x2
      p = p2

    print "finished - all masked"
    indexsorted = numpy.argsort(xminlist)
    xminsorted = []
    pminsorted = []
    for n in indexsorted :
      xminsorted.append( xminlist[n] )
      pminsorted.append( pminlist[n] )
    print numpy.array_str( numpy.array(xminsorted), precision=4 ) 
    gap = []
    for n in range(1,len(xminsorted)) :
      gap.append( xminsorted[n] - xminsorted[n-1] )
    print numpy.array_str( numpy.array(gap), precision=4 ) 
    return [xminsorted, pminsorted]

def processFile( infile ) :

    [x,p,fGHz] = readscan( infile )
    print "fGHz = %.2f; lambda/2 = %.4f mm" % (fGHz,  0.5*299.8/fGHz)

    x1 = 0.
    while (x1 >= 0.) :
      xminsorted,pminsorted = findminima( x, p, fGHz )
      print "enter -1 to get out of loop"
      startstop = raw_input( "startpos stoppos (mm): ") 
      a = startstop.split()
      try :
        x1 = float(a[0])
      except:
        x1 = 0.
      try :
        x2 = float(a[1])
      except :
        x2 = 20.
      print "averaging values from %.3f to %.3f mm" % (x1,x2)
      navg = 0
      avg = 0.
      for i in range(0,len(x)) :
        if (x[i] >= x1) and (x[i] <= x2) :
            avg = avg + p[i]
            navg = navg + 1
      print "%s: n = %d, avg = %.5f" % (infile, navg, avg/navg )
      print " "

    # append values to csv file
      fout = open( "summary2.csv", "a" )
      fout.write("%s,%.4f,%.4f,%d,%.5f\n" % (infile, x1, x2, navg, avg/navg))
      fout.close()

def readscan( infile ) :
    xin = []
    pin = []
    try :
      n = infile.find("_")
      fGHz = float( infile[0:n] )
    except :
      fGHz = float( raw_input( "fGHz: ") )
    fin = open( infile, "r" )
    for line in fin :
      if line.startswith("#") :
        print line
      else :
        a = line.split()
        xin.append( float(a[1]) )
        pin.append( float(a[2]) )
    return numpy.array(xin), numpy.array(pin), fGHz
    
def process( scanList=["100_out","100_posA","100_posB"], plotName="100" ) :
    pyplot.ioff()
    pp = PdfPages("%s.pdf" % plotName )
    npanels = 3
    fig = pyplot.figure( figsize=(8,11) )

    npanel = 0
    for scan in scanList: 
      npanel = npanel + 1
      [x,p,fGHz] = readscan( scan )
      s = pyplot.subplot(npanels,1,npanel)
      s.grid( True, linewidth=0.05, color="0.01" )   # color=0.1 is a light gray

      if npanel == 1 :
        xmin = x.min() - 0.04*(x.max()-x.min()) 
        xmax = x.max() + 0.04*(x.max()-x.min()) 
      pmin = p.min() - 0.04*(p.max()-p.min()) 
      pmax = p.max() + 0.06*(p.max()-p.min()) 

      s.axis( [xmin, xmax, pmin, pmax] )
      s.plot( x, p, color="r", linewidth=0.5 )
      xminsorted,pminsorted = findminima( x, p, fGHz )
      s.plot( xminsorted, pminsorted, "bo" )
      s.text(.03, .92,  "%s" % scan, horizontalalignment='left', transform=s.transAxes, fontsize=12 )

    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()

def calib( infile, outfile ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  fout.write("# calibration from %s\n" % infile)
  first = True
  for line in fin :
    a = line.split()
    if a[0] == "#" :
      dB = float( a[1] )
      if first :
        dBref = dB  
      v = []
    else :
      v.append( float( a[0] ))
      if len(v) >= 2 :
        vavg = (v[0] + v[1])/2.
        trueTrans = pow(10., -(dB-dBref)/10.)
        if first :
          vref = vavg 
          first = False 
        measRatio = vavg/vref
        print "%.2f   %.5f   %.5f" % (dB,trueTrans,measRatio)
        fout.write( "%.2f   %.5f   %.5f\n" % (dB,trueTrans,measRatio) )
  fout.close()

def doit() :
  process( ["70_out1", "70_posA", "70_posB"], "70a" )
  process( ["70_out2", "70_posA", "70_posB"], "70b" )
  process( ["80_out", "80_posA", "80_posB"], "80" )
  process( ["90_out", "90_posA", "90_posB"], "90" )
  process( ["100_out", "100_posA", "100_posB"], "100" )
  process( ["110_out", "110_posA", "110_posB"], "110a" )
  process( ["110_outrepeat", "110_posA", "110_posB"], "110b" )

def doit2() :
  calib( "cal_70", "calib_70" )
  calib( "cal_80", "calib_80" )
  calib( "cal_90", "calib_90" )
  calib( "cal_100", "calib_100" )
  calib( "cal_110", "calib_110" )

def doit3() :
  processFile( "70_out1" )
  processFile( "70_out2" )
  processFile( "70_posA" )
  processFile( "70_posB" )
  processFile( "80_out" )
  processFile( "80_posA" )
  processFile( "80_posB" )
  processFile( "90_out" )
  processFile( "90_posA" )
  processFile( "90_posB" )
  processFile( "100_out" )
  processFile( "100_posA" )
  processFile( "100_posB" )
  processFile( "110_outrepeat" )
  processFile( "110_out" )
  processFile( "110_posA" )
  processFile( "110_posB" )
