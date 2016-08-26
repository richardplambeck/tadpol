# leo.py
# routines for calculating and plotting galaxy rotation curves
# plotpv2 adapted from ori3.plotpv 

import math
import time
import cmath
import numpy
import sys
import pickle
import subprocess
import shlex
import string
import random
import matplotlib
import os
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import signal

# ========================== plot position velocity diagram ========================== #

def extractParameter( instring, param ) :
    '''extract parameter "param" from imlist output, passed to this routine as instring'''
    n = instring.find( param )
    if n > 0 :
      a = instring[n:].split()
      return float( a[2] )
    else :
      return None
  
# ----------------------------------------------------------------------------------------- #
# process pv diagram written out by velplot as a miriad map
# axis1 is VELO; axis2 is position; axis3 is intensity
# 4/20/16 added vmin,vmax,pmin,pmax as a quick way of selecting piece of image that I want;
#   hopefully defaults are so large that they won't trip me up in cases where I don't care

def readpv( imageFile, vmin=-1.e8, vmax=1.e8, pmin=-1.e8, pmax=1.e8 ) :
    '''read pos-vel map in Miriad format created by velplot, convert to numpy array'''

  # begin by figuring out axes
    p = subprocess.Popen( ( shlex.split('imlist in=%s units=absolute' % imageFile ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    vstart = extractParameter( result, "crval1" )
    pstart = extractParameter( result, "crval2" )
    vstep = extractParameter( result, "cdelt1" )
    pstep = extractParameter( result, "cdelt2" )
    print "vstart = %.2f, vstep = %.2f, pstart = %.3f, pstep = %.3f" % (vstart,vstep,pstart,pstep) 

  # now dump the data and assign to array
    p = subprocess.Popen( ( shlex.split("imtab in=%s log=imtablog format=(3F12.5)" % imageFile ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)  
    result = p.communicate()[0]
    print result

  # open the logfile, fill out lists, then convert to array
    vel = []
    pos = []
    flx = []
    fin = open("imtablog", "r")
    for line in fin :
      a = line.split()
      v = float(a[0]) + vstart
      p = float(a[1]) + pstart
      if (v >= vmin) and (v <= vmax) and (p >= pmin) and (p <= pmax) :
        vel.append( v )   
        pos.append( p )
        flx.append( float(a[2]) )
    fin.close()
    return [ numpy.array(vel), numpy.array(pos), numpy.array(flx) ]

# this is a special purpose routine to search selected directories for the pv.mp file
#   for spectral line with frequency linefreq
# the prefix to the file can contain any character string, but file must end in
#   "xxx.xx.pv.mp", where xxx.xx is linefreq to 2 decimal places

def findpvFile( linefreq, searchType=".pv.mp" ) :
    dirList = ["/big_scr6/plambeck/350GHz/miriad", "/big_scr6/plambeck/650GHz/miriad"]
      # list of directories that will be searched
    searchstring = "%.2f" % linefreq + searchType
      # filename must include linefreq in GHz, rounded to 2 decimal accuracy
    for directory in dirList :
      for filename in os.listdir( directory ) :
        if searchstring in filename :
          return directory + "/" + filename
    print "couldn't find filename with searchstring %s" % searchstring
    return None

# pv plot for LeoCO
def plotpv2( imageFile, cutoff ) :
    pyplot.ioff()
    pp = PdfPages("spectrum.pdf")
    fig = pyplot.figure( figsize=(8,11) )

  # plot unmasked file at top
    ax1 = fig.add_axes( [.1,.6,.7,.3])      # main plot
    ax2 = fig.add_axes([.85,.6,.015,.3])    # color wedge
    ax2.tick_params( labelsize=10 )
    ax1.tick_params( labelsize=10 )
    [vel, pos, flx] = readpv( imageFile )
       # velocity changes fastest; one pos row per velocity scan
    print vel
    print pos
    print flx
    nv = len( numpy.unique( vel ) )
    np = len( numpy.unique( pos ) )
    print "nv = %d, np = %d" % (nv,np)
    print pos[:-1:nv]

  # be careful here - labeling can very easily be reversed
  # first shape array as (np,nv) - 2nd index varies fastest; effectively, the flx vector is
  #   broken into np vectors of length nv, which is what we want 
  # then transpose the array so that it consists of nv rows of np elements; this is more
  #   suitable for plotting because there are more positions than velocities
    flx2 = numpy.transpose(numpy.reshape(flx, (np,nv) ))

  # plotting with 'origin='upper' puts element 0,0 in the upper left corner; 
  # extents must be labeled pos[0],pos[-1],vel[-1],vel[0] to match
    imgplot = ax1.imshow(flx2, origin='upper', aspect='auto', \
        extent=[pos[0],pos[-1],vel[-1],vel[0]] )

    ax1.grid( True, linewidth=1, color="white")   # color=0.1 is a light gray
    pyplot.colorbar( imgplot, cax=ax2  )

  # plot masked file at bottom
    ax1 = fig.add_axes( [.1,.25,.7,.3])      # main plot
    ax2 = fig.add_axes([.85,.25,.015,.3])    # color wedge
    ax2.tick_params( labelsize=10 )
    ax1.tick_params( labelsize=10 )
    palette = pyplot.get_cmap("jet")
    palette.set_bad('gray',1.0)
    flx3 = numpy.ma.masked_where( flx2 < cutoff, flx2) 
    imgplot = ax1.imshow(flx3, origin='upper', aspect='auto', \
        extent=[pos[0],pos[-1],vel[-1],vel[0]], cmap=palette )
    ax1.grid( True, linewidth=1, color="white")   # color=0.1 is a light gray
    pyplot.colorbar( imgplot, cax=ax2  )

  # use weighted average to get rotation curve
    print vel[0:nv]
    print flx3[:,[10]]   # expect 41 outputs
    vcurve = numpy.zeros( np, dtype=float)
    #print numpy.ma.sum( flx3[0:nv][40] )
    #print numpy.ma.dot( vel[0:nv], flx2[0:nv][40] )/numpy.ma.sum( flx2[0:nv][40] )
    for npp in range(0,np) :
       vcurve[npp] = numpy.ma.dot( vel[0:nv], flx3[:,[npp]] )/numpy.ma.sum( flx3[:,[npp]] )
    print vcurve
    ax1.plot( pos[:-1:nv], vcurve, 'b-')
    ax1.axis( [pos[0], pos[-1], vel[-1], vel[0]] )

    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show() 
    #print vel
    #print pos
    #print flx 
      

