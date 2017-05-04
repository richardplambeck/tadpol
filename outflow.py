# def outflow.py
# I am starting new collection of functions to model orion outflow

import math
import time
import numpy
import subprocess
import shlex
import string
import sys
import os
import matplotlib
matplotlib.use('GTKAgg')
from matplotlib import pyplot
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.widgets import Cursor
from matplotlib.mlab import PCA


startx = 0.
starty = 0.

# onclick is called by fig.canvas.mpl_connect when a mouse button is clicked
# I am using it to identify streamers or filaments in orion maps (mostly CO)
# I can't seem to get onclick to communicate via global variables (defined outside
#    any routine), so I am writing data to file "mouseclicks.dat" each time the 
#    the mouse is clicked
# left mouse button (event.button = 1) is used to identify beginning of streamer
# right mouse button (event.button = 3) is used to identify the end 
#
def onclick(event) :
    fout = open("mouseclicks.dat","a")
    fout.write(" %d  %.3f  %.3f\n" % (event.button, event.xdata, event.ydata))
    fout.close()
      
# returns starting velocity, velocity step, number of velocities in an image cube
def getvelocitydata( imagefile ) : 
    p= subprocess.Popen( ( shlex.split('imlist in=%s' % imagefile) ), \
        stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    lines = result.split("\n")
    for line in lines :
      if len(line) > 1 :
        a = line.split()
        n = string.find( line, "naxis3  :" )
        if n >= 0 :
          nchans = int( line[n+9:].split()[0] )
        n = string.find( line, "crval3  :" )
        if n >= 0 :
          v1 = float( line[n+9:].split()[0] )
        n = string.find( line, "cdelt3  :" )
        if n >= 0 :
          dv = float( line[n+9:].split()[0] )
    #print "nchans = %d; v1 = %.3f km/sec; dv = %.3f km/sec" % (nchans,v1,dv)
    return nchans, v1, dv
    
# this is copied from marsPA.py; dumps data from Miriad to logfile, then read log file to create numpy arrays
def readArray( imageFile, region="arcsec,box(19,-30,-50,40)", image=32 ) :
    p = subprocess.Popen( ( shlex.split('imtab in=%s region=%s(%d) log=imtablog' % ( imageFile, region, image ) ) ), \
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
    #p.terminate()
    return [ numpy.array(x), numpy.array(y), numpy.array(z) ]

def plotArray( x,y,z,velocity ) :
    #pyplot.ion()
    fig = pyplot.figure()
    p = fig.add_axes( [.08,.05,.8,.8])
    p2 = fig.add_axes([.91,.05,.015,.35])
    xu = numpy.unique(x)    # sorted unique elements in x array
    nx = len(xu)
    yu = numpy.unique(y)    # sorted unique elements in y array
    ny = len(yu)
    imgplot = p.imshow( numpy.reshape(z,(ny,nx)), origin='lower', \
      extent=[xu[0],xu[-1],yu[0],yu[-1]] )
    #p.set_xlim( xu[0],xu[-1] )
    #p.set_ylim( yu[0],yu[-1] )
    pyplot.colorbar( imgplot, cax=p2  )
    cursor = Cursor( p, useblit=True, color='red', linewidth=1 )
    cid = fig.canvas.mpl_connect('button_press_event', onclick )  
    parseclicks( p, "mouseclicks.dat", velocity )
    pyplot.show()

# display each velocity channel, use mouse to identify streamers (left
#	button for beginning of streamer, right button for the end); when
#	finished with velocity channel, close the image on screen - next
#   image will be plotted
# enter velocity in the file "mouseclicks.dat" each time a new chan is plotted
def markstreamers( imagefile ) :
    nchans, v1, dv = getvelocitydata( imagefile ) 
    print "step through %d images" % nchans
    #fout = open("mouseclicks.dat","a")
    #fout.write(" image: %s\n" % imagefile )
    #fout.close()
    for n in range(0,nchans) :
      x,y,z = readArray( imagefile, image=n+1 )
      velocity = v1+n*dv
      print "## VELOCITY = %.2f ##" % velocity
      fout = open("mouseclicks.dat","a")
      fout.write(" velocity: %.3f\n" % velocity )
      fout.close()
      plotArray( x,y,z,velocity )

# read through mouseclicks.dat, draw line segments for each streamer
def parseclicks( p, infile, v ) :
    x = [-1000.,-1000.]
    y = [-1000.,-1000.]
    fin = open(infile, "r")
    for line in fin :
      n = string.find( line, "velocity:" )
      if n>= 0 :
        velocity = float( line[n+9:] )
        if velocity == v :
          useSegments = True
        else :
          useSegments = False 
      else :
        a = line.split()
        if int(a[0]) == 1 :
          x[0] = float(a[1])
          y[0] = float(a[2])
        elif (int(a[0]) == 3) and (x[0] != -1000.) and useSegments :
          x[1] = float(a[1])
          y[1] = float(a[2])
          p.plot(x,y,color='black')
          x[0] = -1000.

#imagefileList = ["SO_219.95.cm","SO2_216.64.cm","SiO_217.10.cm"] 
#imagefileList = ["CO_230.54.cm", "SiO_217.10.cm"] 
#imagefileList = ["SiS_217.82.cm","SiO_217.10.cm"] 
#imagefileList = ["SiO_217.10.cm", "SO_219.95.cm", "HC3N_218.32.cm"]
#imagefileList = ["SO_219.95.cm","SO2_216.64.cm","CH3OH_231.28.cm","H2S_216.71.cm","SiS_217.82.cm","HC3N_218.32.cm"] 
imagefileList = ["SO_219.95.avg.cm","SO2_216.64.avg.cm","CH3OH_231.28.avg.cm","H2S_216.71.avg.cm", \
                  "SiS_217.82.avg.cm","HC3N_218.32.avg.cm","CO_230.54.avg.cm", "SiO_217.10.avg.cm"]
#imagefileList = ["CO_230.54.avg.cm","SiO_217.10.avg.cm"]
#imagefileList = ["SiO_217.10.avg.cm","SiS_217.82.avg.cm"]
def mypca( imagefileList=imagefileList, n=20 ) :
    first = True
    nimages = 0
    for imagefile in imagefileList :
      nimages = nimages+1
      print imagefile
      #nchans, v1, dv = getvelocitydata( imagefile ) 
      #for n in range(0,nchans) 
      x,y,z = readArray( imagefile, image=n+1 )
      if first :
        d = z
        first = False
      else :
        d = numpy.concatenate( (d,z) )
    nrows = len(d)/nimages
    d = numpy.swapaxes( numpy.reshape(d,(-1,nrows)), 0, 1 )
    myPCA = PCA(d)
    print "fracs: ", myPCA.fracs
    z = myPCA.Y[:,0]
    fig = pyplot.figure()
    p = pyplot.subplot( 1,2,1)
    xu = numpy.unique(x)    # sorted unique elements in x array
    nx = len(xu)
    yu = numpy.unique(y)    # sorted unique elements in y array
    ny = len(yu)
    imgplot = p.imshow( numpy.reshape(z,(ny,nx)), origin='lower', \
      extent=[xu[0],xu[-1],yu[0],yu[-1]] )
    p = pyplot.subplot( 1,2,2)
    z = myPCA.Y[:,1]
    imgplot = p.imshow( numpy.reshape(z,(ny,nx)), origin='lower', \
      extent=[xu[0],xu[-1],yu[0],yu[-1]] )
    pyplot.show()

def panels( imagefileList=imagefileList, region='arcsec,box(19,-30,-50,40)', image=0 ) :
    fig = pyplot.figure()
    np = 0
    for imagefile in imagefileList :
      np = np+1
      p = pyplot.subplot( 2,4,np )
      x,y,z = readArray( imagefile, region=region, image=image+1 )
      xu = numpy.unique(x)    # sorted unique elements in x array
      nx = len(xu)
      yu = numpy.unique(y)    # sorted unique elements in y array
      ny = len(yu)
      print xu[0],xu[-1]
      imgplot = p.imshow( numpy.reshape(z,(ny,nx)), origin='lower', \
        extent=[xu[0],xu[-1],yu[0],yu[-1]] )
      p.set_xlim( xu[0],xu[-1] )
      p.set_ylim( yu[0],yu[-1] ) 
      p.set_title( imagefile, fontsize=10 )
      p.plot( 0., 0., "+", color="white", markersize=4., markeredgewidth=2. )
    pyplot.show()
  
