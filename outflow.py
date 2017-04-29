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


def onclick(event) :
    fout = open("mouseclicks.dat","a")
    fout.write(" %d  %.3f  %.3f\n" % (event.button, event.xdata, event.ydata))
    fout.close()

# this is copied from marsPA.py; dumps data from Miriad to logfile, then read log file to create numpy arrays
def readArray( imageFile, region="arcsec,box(19,-30,-50,40)(32)" ) :
  p = subprocess.Popen( ( shlex.split('imtab in=%s region=%s log=imtablog' % ( imageFile, region ) ) ), \
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
  #return [ numpy.array(x), numpy.array(y), numpy.array(z) ]

#def plotArray( x,y,z ) :
  fig = pyplot.figure()
  p = fig.add_axes( [.08,.05,.8,.8])
  p2 = fig.add_axes([.91,.05,.015,.35])
  xu = numpy.unique(x)    # sorted unique elements in x array
  nx = len(xu)
  yu = numpy.unique(y)    # sorted unique elements in y array
  ny = len(yu)
  print nx,ny
  imgplot = p.imshow( numpy.reshape(z,(ny,nx)), origin='lower', \
    extent=[xu[0],xu[-1],yu[0],yu[-1]] )
  #p.set_xlim( xu[0],xu[-1] )
  #p.set_ylim( yu[0],yu[-1] )
  pyplot.colorbar( imgplot, cax=p2  )
  cursor = Cursor( p, useblit=True, color='red', linewidth=2 )
  cid = fig.canvas.mpl_connect('button_press_event', onclick)  
  pyplot.show()

