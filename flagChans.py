# flagChans.py
# this is designed to find, and flag, spectral lines in OriALMA data
# 
# procedure:
#  - create text file with sectrum with uvspec:
#      uvspec vis=xx axis=chan,amp options=avall,nobase,ampscalar interval=10000 \
#         select=... stokes=I log=blah
#      or, create vector average if desired
#  - edit "doit" to plot these files
#  - examine files on screen, edit flag file showing what channels to flag;
#      NOTE: flags are in miriad numbering convention (chans 1-nch), not (0 to nchan-1)
#  - process 500 channels at a time; find 10th percentile in amplitude,
#      flag everything above this
#  - plot results

import numpy
import subprocess
import shlex
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot

# readspec reads text file created by uvspec (usually with options=ampscalar,avall,nobase)
def readspec( infile, flagtable=None ) :
    chn = numpy.arange(0,1921,1) 
    amp = numpy.full( 1920, numpy.nan )
    fin = open( infile, "r" )
    for line in fin:
      a = line.split()
      n = int( round( float(a[0]) ) )       # round accepts and returns a float
      amp[n-1] = float(a[1])
    fin.close()
    # return numpy.array(amp)
    print amp
    fitbaseline( numpy.array(amp), outfile='junk', nchPerSegment=192, flagtable=flagtable ) 
   
# fitbaseline splits amplitude array into chunks, finds minimum value in each chuck,
def fitbaseline( amp, outfile='junk', nchPerSegment=192, flagtable=None ) :
    nch = len(amp)    # number of channels in array
    ndiv = int( round(nch/nchPerSegment ) )
    print "... spectrum has %d channels" % nch
    print "... splitting spectrum into %d segments, %d channels per segment" % (ndiv, nch/ndiv)
    segments = numpy.split( amp, ndiv )
      # split amp into ndiv sections
    nSegmin = numpy.nanargmin(segments, axis=1)
      # return array containing channel number of min within each segment
    nBase = []
    yBase = []
    for n in range(0,ndiv) :
      nBase.append( n*(nch/ndiv) + nSegmin[n] )
        # channel numbers of minima in amp array
      yBase.append( amp[ nBase[-1] ])
    xBase = numpy.array( nBase, dtype=float )
    xInterp = numpy.zeros( len(amp) )
    for n in range(0,len(amp)) :
      xInterp[n] = float(n)
    yInterp = numpy.interp( xInterp, xBase, yBase )
    fout = open( outfile, "w" )
    for n in range(0,len(xInterp)) :
      fout.write("%5d %8.5f\n" % (n,yInterp[n]) )
    fout.close()
    ratio = amp/yInterp
    r = numpy.ma.masked_greater( ratio, 1.05 )
    legity = numpy.ma.masked_where( numpy.ma.getmask(r), amp)
    legitx = numpy.ma.masked_where( numpy.ma.getmask(r), xInterp)
    pyplot.ioff()
    fig = pyplot.subplot( 1,1,1 )
    print "opening plot"
    fig.axis( [xInterp[0],xInterp[-1],numpy.nanmin(amp), numpy.nanmax(amp)] )
    fig.plot( xInterp, amp, "-",color='blue')
    fig.plot( xInterp, yInterp, "-", color="red")
    fig.plot( legitx, legity, "-", color='cyan')
    pyplot.minorticks_on()
    fig.grid( which='both' )
    #fig.grid( )

  # indicate channels that were flagged bad for continuum map
    if flagtable :
      y0 = numpy.zeros( len(xInterp) )
      flaggedBad = numpy.zeros( len(xInterp), dtype=bool )
      fin = open( flagtable, "r" )
      for line in fin :
        if not line.startswith("#") :
          a=line.split(",")
          nbstart = int(a[0]) - 1
          nbstop = int(a[1]) - 1
          for n in range(0, len(xInterp)) :
            if (xInterp[n] >= nbstart) and (xInterp[n] <= nbstop) :
              flaggedBad[n] = True
      fin.close()
      fig.fill_between( xInterp, amp, 0., where=flaggedBad, color='red', alpha=0.2 )
    print "show plot"
    pyplot.show()

# working from list of bad chans, print out uvflag commands to be copied into master script
def flagCmds( maskFile, spw="$WIN.uv" ) :
  fin = open( maskFile, "r" )
  for line in fin :
    if not line.startswith("#") :
      a = line.split(',')
      nstart = int(a[0]) 
      nstop = int(a[1]) 
      print "    uvflag vis=%s flagval=flag line=chan,%d,%d,1" % (spw,nstop-nstart+1,nstart)

# I could also do this directly from python, but maybe this is not a good idea
#  p= subprocess.Popen( ( shlex.split('gplist vis=%s options=all' % infile) ), \
#     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
#  result = p.communicate()[0]


def doit() :  
    # readspec( 'spw0.gt200.ampvector.txt', flagtable="flags_a0" )
    # readspec( 'spw1.gt200.ampvector.txt', flagtable="flags_a1" )
    # readspec( 'spw2.gt200.ampvector.txt', flagtable="flags_a2" )
    # readspec( 'spw3.gt200.ampvector.txt', flagtable="flags_a3" )
    readspec( 'spw0.gt200.ampscalar.txt', flagtable="flags_a0" )
    readspec( 'spw1.gt200.ampscalar.txt', flagtable="flags_a1" )
    readspec( 'spw2.gt200.ampscalar.txt', flagtable="flags_a2" )
    readspec( 'spw3.gt200.ampscalar.txt', flagtable="flags_a3" )
    #readspec( 'BN.spw0.gt200.ampscalar.txt', flagtable='BNflags_b0' ) 
    #readspec( 'BN.spw1.gt200.ampscalar.txt', flagtable='BNflags_b1' )
    #readspec( 'BN.spw2.gt200.ampscalar.txt', flagtable='BNflags_b2' )
    #readspec( 'BN.spw3.gt200.ampscalar.txt', flagtable='BNflags_b3' )

