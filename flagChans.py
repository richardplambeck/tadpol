# flagChans.py
# this is designed to find, and flag, spectral lines in OriALMA data
# 
# procedure:
#  - read in spectrum
#  - process 500 channels at a time; find 10th percentile in amplitude,
#      flag everything above this
#  - plot results

import numpy
import subprocess
import shlex
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot

def readspec( infile, flagtable=None ) :
  amp = []
  fin = open( infile, "r" )
  for line in fin:
    a = line.split()
    amp.append(float(a[1]))
  fin.close()
  # return numpy.array(amp)
  fitbaseline( numpy.array(amp), outfile='junk', nchPerSegment=192, flagtable=flagtable ) 
   
def fitbaseline( amp, outfile='junk', nchPerSegment=192, flagtable=None ) :
  ndiv = int( round(len(amp)/nchPerSegment ) )
  print "... splitting spectrum into %d segments, %d channels per segment" % (ndiv, len(amp)/ndiv)
  segments = numpy.split( amp, ndiv )
	# split amp into ndiv sections
  nSegmin = numpy.argmin(segments, axis=1)
    # return array containing channel number of min within each segment
  nBase = []
  yBase = []
  for n in range(0,ndiv) :
    nBase.append( n*(3840/ndiv) + nSegmin[n] )
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
  pyplot.ion()
  fig = pyplot.subplot( 1,1,1 )
  print "opening plot"
  fig.axis( [xInterp[0],xInterp[-1],amp.min(),amp.max()] )
  fig.plot( xInterp, amp, "-",color='blue')
  fig.plot( xInterp, yInterp, "-", color="red")
  fig.plot( legitx, legity, "-", color='cyan')
  #fig.grid( b=True, which='minor' )
  fig.grid( )

# indicate channels that were flagged bad for continuum map
  if flagtable :
    y0 = numpy.zeros( len(xInterp) )
    flaggedBad = numpy.zeros( len(xInterp), dtype=bool )
    fin = open( flagtable, "r" )
    for line in fin :
      if not line.startswith("#") :
        a=line.split(",")
        nbstart = int(a[0])
        nbstop = int(a[1])
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


  
