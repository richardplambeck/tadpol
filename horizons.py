# goal is to generate table of Jupiter moon offsets and fluxes from JPL horizons emails

import subprocess
import shlex
import math
import numpy
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

clight = 2.998e10  # cm/sec
kB = 1.38e-16      # erg/K
Jy = 1.e-23        # erg/(sec-cm^2-Hz)
h = 6.626e-27

def SJy( diamSec, TK, freqGHz=150. ) :
  dOmega = math.pi * pow( diamSec*math.pi/(2.*180.*60.*60.), 2.)    # steradians
  SJy =  2.* h * pow( freqGHz*1.e9, 3.) / pow( clight, 2.)  * 1./(math.exp( h * freqGHz*1.e9/( kB * TK)) - 1.) * dOmega / Jy
  return SJy


def parseFile( infile, outfile ):
    fin = open( infile, "r" )
    fout = open( outfile, "w" )
    parseLine = False
    for line in fin :
      if line.startswith( "$$SOE" ) :
        parseLine = True
      elif line.startswith( "$$EOE" ) :
        parseLine = False
      else : 
        if parseLine :
          a = line.split()
          nf = 0
          calDate = a[nf] 
          calTime = a[nf+1] 
          JD = float( a[nf+2] ) 

        # sometimes date is followed by non-blank solar or lunar flags
        # increment nf to skip over these if they are present
          try: 
            az = float( a[nf+3] )   # fails if not a float
          except :
            nf = nf+1
            az = float( a[nf+3] )
          el = float( a[nf+4] ) 
        
        # separation is followed by '/c' where c is a code
          slashpos = a[nf+5].find('/')
          sep = float( a[nf+5][0:slashpos] ) 
          diam = float( a[nf+6] ) 
          fout.write("%12s %5s %12.3f %10.4f %8.4f %8.3f %6.3f\n" % \
            (calDate, calTime[0:5], JD, az, el, sep, diam))
    fin.close()
    fout.close()

def readfile( infile ) :
    calDate = []
    calTime = []
    JD = []
    az = []
    el = []
    sep = [] 
    diam = []
    fin = open( infile, "r" )
    for line in fin :
      a = line.split()
      calDate.append( a[0] )
      calTime.append( a[1] )
      JD.append( a[2] )
      az.append( float( a[3] ) )
      el.append( float( a[4] ) )
      sep.append( float( a[5] ) )
      diam.append( float( a[6] ) )
    fin.close()
    return [calDate,calTime,JD,numpy.array(az),numpy.array(el),numpy.array(sep),numpy.array(diam)]

def makeSepTable( jupfile, moonfile, Tb, miriadfile ) :
    jDate,jTime,jJD,jaz,jel,jsep,jdiam = readfile( jupfile )
    mDate,mTime,mJD,maz,mel,msep,mdiam = readfile( moonfile )
    delaz = 60. * (maz - jaz) * numpy.cos(numpy.radians(jel))
    delel = 60. * (mel-jel )
    sepcalc = 60.* numpy.sqrt( pow(delaz,2) + pow(delel,2))
    fout = open( miriadfile, "a" )
    for n in range(0, len(jaz)) :
      print mJD[n], delaz[n], delel[n], sepcalc[n], msep[n]
      fout.write("imgen in=test.mp out=test1.mp factor=1 object=gaussian spar=%.2f,%.2f,%.2f,180.,180.,0. options=totflux\n" %
                        ( SJy( mdiam[n], Tb, 150.), 60.*delaz[n], 60.*delel[n] ) )
      fout.write("rm -fr test.mp\n")
      fout.write("mv test1.mp test.mp\n")
    fout.close()
    return delaz, delel    

def makePlot() :
    pyplot.ioff()
    pp = PdfPages("moons.pdf")
    fig = pyplot.figure()
    ax = pyplot.subplot(1,1,1)
    #ax.axis( [-20.,20.,-20.,20.] )
    ax.axis( [-15.,15.,-15.,15.] )
    ax.set_aspect('equal')
    ax.plot( [0.], [0.], '+')
    ax.grid( True, linewidth=0.05, color='0.1' )
    #for moon,color in zip( ["Callisto"],["blue"] ) :
    makeSepTable( "Jupiter.azel", "Jupiter.azel", 165., "miriadcmd" )
       # this is just to put jupiter on the Miriad map
    for moon,Tb,color in zip( ["Callisto","Europa","Ganymede","Io"],[112.,94.,92.,98.],["blue","green","orange","red"] ) :
       delaz,delel = makeSepTable( "Jupiter.azel", "%s.azel" % moon, Tb, "miriadcmd" )
       ax.plot( delaz, delel, marker="o", linestyle='None', color=color, label=moon)
    ax.legend( loc=0, prop={'size':10}, numpoints=1 )
    ax.set_xlabel('delta az * cos(el) (arcmin)')
    ax.set_ylabel('delta el (arcmin)')
    pyplot.title('Jupiter moon positions during season 1 observations', fontdict={'fontsize':12} ) 
    pyplot.savefig( pp, format='pdf')
    pyplot.show()
    pp.close()   

def formatQuery( infile='1stSeason_JupiterObservations.txt', outfile='tlist' ) :
    fout = open( outfile, "w" )

  # load the basic query stuff in the outfile
    fout.write("!$$SOF\n")
    fout.write("COMMAND = '599' '501' '502' '503' '504'\n")
        # object codes: 599=Jupiter, 501=Io, 502=Europa, 503=Ganymede, 504=Callisto
    fout.write("CENTER = '-7'\n")
        # topocentric results for ALMA array
    fout.write("OBJ_DATA = 'NO'\n")
    fout.write("TABLE_TYPE = 'OBS'\n")
    fout.write("QUANTITIES = '4,12,13'\n")
        # 4 = apparent az and el, 12 = satellite offset (arcsec), 13 = target angular diam
    fout.write("CAL_FORMAT = BOTH\n")
    fout.write("EXTRA_PREC = 'YES'\n")
        # return az el (degrees) to 4 significant figures

  # now generate TLIST
    fout.write("TLIST =\n")
    fin = open( infile, "r" )
    for line in fin:
      a = line.split()
      tstart = float( a[0] ) 
      tstop = float( a[1] ) 
      tmed = (tstart + tstop)/2. #+ 2400000.5
      print "   DELTA = %.3f" % (tstop - tstart)
      
      fout.write(" '%.5f'\n" % tmed )
    fin.close()
    fout.write("\n! START_TIME='2014-JAN-31 16:00'\n")
    fout.write("!$$EOF")
    fout.close()   
 
def doit() :
    parseFile( 'Jupiter.dat', 'Jupiter.azel' )
    parseFile( 'Io.dat', 'Io.azel' )
    parseFile( 'Europa.dat', 'Europa.azel' )
    parseFile( 'Ganymede.dat', 'Ganymede.azel' )
    parseFile( 'Callisto.dat', 'Callisto.azel' )
    makePlot()   

# adapted from marsPA.py
def readArray( imageFile ) :
  # dump selected region of image to a logfile
  p = subprocess.Popen( ( shlex.split('imtab in=%s log=imtablog' % imageFile ) ), \
     stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
  result = p.communicate()[0]
  print result
  # open the logfile, fill out lists
  x = []
  y = []
  z = []
  fin = open("imtablog", "r")
  for line in fin :
    a = line.split()
    x.append( -1.*float(a[0]) )   # reverse RA to get arcsec on sky
    y.append( float(a[1]) )
    z.append( float(a[2]) )
  fin.close()
  return [ numpy.array(x), numpy.array(y), numpy.array(z) ]

# reads in Miriad model
def beamplot( imageFile ='test.mp' ) :
    pyplot.ioff()
    pp = PdfPages("map.pdf")
    fig = pyplot.figure( figsize=(8,8) )

  # plot unmasked file at top
    ax1 = fig.add_axes( [.1,.1,.7,.7])      # main plot
    ax2 = fig.add_axes([.85,.1,.03,.7])    # color wedge
    ax2.tick_params( labelsize=12 )
    ax1.tick_params( labelsize=12 )
    [x,y,z] = readArray( imageFile )
    x = x/60.
    y = y/60.
    zmax = z.max()
    print "zmax = %.2e" % zmax
    z1 = 10.*numpy.log10(z/zmax)
    z2 = numpy.ma.masked_where( z1 < -50., z1 )
    ax1.axis( [x[0], x[-1], y[0], y[-1]] )
    nx = len( numpy.unique( x ) )
    ny = len( numpy.unique( y ) )
    print "nx = %d, ny = %d" % (nx,ny)
    z3 = numpy.reshape( z2, (nx,ny) )
    
    #pyplot.pcolor( x, y, z3, norm=LogNorm(vmin=1.e-5,vmax=1.))
    #pyplot.colorbar()

    imgplot = ax1.imshow(z3, origin='upper', aspect='auto',\
       extent=[x[-1],x[0],y[-1],y[0]] )
    pyplot.title('dB', fontdict={'fontsize':12} ) 
    ax1.set_xlabel('delta az*cos(el) (arcmin)')
    ax1.set_ylabel('delta el (arcmin)')

    pyplot.colorbar( imgplot, cax=ax2  )

    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()
      

