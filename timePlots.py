# timePlots.py
#
# this is a collection of routines that use the matplotlib.dates module
#   to make plots with time as the x-axis

import math
import sys
import time
import numpy
import subprocess
import shlex
import string
import leakSolve
import paPlot
import pickle
import datetime
import scipy.optimize
import dateutil
import matplotlib.pyplot as pyplot
import matplotlib.dates as mpdates
from matplotlib.backends.backend_pdf import PdfPages

     
def getdata(infile,col) :
  t = []
  s = []
  fin = open(infile,"r")
  for line in fin:
    if not line.startswith("#"):
      a = line.split()
      #t.append(mpdates.date2num(time.strptime(a[1][0:8],"%H:%M:%S")))
      dt = a[0]+":"+a[1][0:8]
      t.append(datetime.datetime.strptime(dt,"%Y-%b-%d:%H:%M:%S"))
      s.append(float(a[col]))
      print dt, a[col] 
  fin.close()
  return[t,s]

def getdata2(infile,col) :
  t = []
  s = []
  fin = open(infile,"r")
  for line in fin:
    if not line.startswith("#"):
      a = line.split()
      t.append(datetime.datetime.strptime( a[0][0:8], "%H:%M:%S") )
      s.append(float(a[col]))
      print a[0],a[col],t,s
  fin.close()
  return[t,s]


def getcsvdata(infile,col) :
  s = []
  fin = open(infile,"r")
  for line in fin:
    if not line.startswith("#"):
      a = line.split(',')
      s.append( a[col] )
      print a[col]
  fin.close()
  return s
  
# convert from LIST of STRINGS containing dates (e.g., "11/9/2011", "11/11/2011" ...] to list of datetime objects
def convert( datestrlist ) :
  out = []
  for dt in datestrlist :
    out.append( dateutil.parser.parse( dt ) )
  return out

def plot( t, s) :
  pyplot.clf()
  pyplot.ion()
  fig = pyplot.subplot(1,1,1)
  tt = mpdates.date2num(t)
  pyplot.ylim(3,42)
  fig.plot_date( tt, s, "-" )
  fig.xaxis.set_major_locator(mpdates.MinuteLocator(interval=10))
  fig.xaxis.set_major_formatter(mpdates.DateFormatter( "%H:%M" ) )
  fig.grid()
  pyplot.ylabel("C11 subreflector temperature (C)")
  pyplot.xlabel("UT on 22-mar-2013")
  #fig.xaxis.set_major_locator(mpdates.AutoDateLocator())
  #fig.xaxis.set_major_formatter(mpdates.AutoDateFormatter())
  pyplot.draw()

#  fig.xaxis.set_major_formatter(timeFmt)
#  fig.draw()
def plot2( t, s) :
  print t,s
  pyplot.clf()
  pyplot.ion()
  fig = pyplot.subplot(2,1,1)
  tt = mpdates.date2num(t)
  #pyplot.ylim(3,42)
  fig.plot_date( tt, s, "-" )
  #fig.xaxis.set_major_locator(mpdates.DateLocator())
  fig.xaxis.set_major_locator(mpdates.AutoDateLocator())
  fig.xaxis.set_major_formatter(mpdates.DateFormatter( "%b%d" ) )
  fig.grid()
  # - get datetime ranges of vlbi schedules and plot them

  ts1 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-21:05:00:00","%Y-%b-%d:%H:%M:%S"))
  ts2 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-21:17:02:00","%Y-%b-%d:%H:%M:%S"))
  pyplot.axvspan(ts1,ts2,facecolor="m",alpha=0.2)

  ts1 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-22:03:00:00","%Y-%b-%d:%H:%M:%S"))
  ts2 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-22:22:07:00","%Y-%b-%d:%H:%M:%S"))
  pyplot.axvspan(ts1,ts2,facecolor="m",alpha=0.2)

  ts1 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-23:03:20:00","%Y-%b-%d:%H:%M:%S"))
  ts2 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-23:22:07:00","%Y-%b-%d:%H:%M:%S"))
  pyplot.axvspan(ts1,ts2,facecolor="m",alpha=0.2)

  ts1 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-25:12:00:00","%Y-%b-%d:%H:%M:%S"))
  ts2 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-25:16:52:00","%Y-%b-%d:%H:%M:%S"))
  pyplot.axvspan(ts1,ts2,facecolor="m",alpha=0.2)

  ts1 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-26:05:00:00","%Y-%b-%d:%H:%M:%S"))
  ts2 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-26:17:06:00","%Y-%b-%d:%H:%M:%S"))
  pyplot.axvspan(ts1,ts2,facecolor="m",alpha=0.2)

  ts1 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-27:12:00:00","%Y-%b-%d:%H:%M:%S"))
  ts2 = mpdates.date2num(datetime.datetime.strptime("2013-Mar-27:18:15:00","%Y-%b-%d:%H:%M:%S"))
  pyplot.axvspan(ts1,ts2,facecolor="m",alpha=0.2)

  pyplot.ylabel("tau230")
  pyplot.xlabel("UT date 2013")
  #fig.xaxis.set_major_formatter(mpdates.AutoDateFormatter())
  pyplot.draw()

def plot3( fig, t, s, serr, ymin, ymax, label=None ) :
  #pyplot.clf()
  #pyplot.ion()
  #fig = pyplot.subplot(2,1,panel)
  tt = mpdates.date2num(t)
  print tt.min(), tt.max()
  ttrange = tt.max() - tt.min()
  fig.axis( [tt.min() - 0.05*ttrange, tt.max() + .05*ttrange, ymin, ymax], size=3 )
  fig.xaxis_date()
  #print x
  #print s
  #print serr

  if len(serr) == len(s) :
    fig.errorbar( tt, s, yerr=serr, fmt='ro', markersize=3 )
  fig.plot( tt, s, "o" )
   
  ##fig.xaxis.set_major_locator(mpdates.DateLocator())
  fig.xaxis.set_major_locator(mpdates.AutoDateLocator())
  fig.xaxis.set_major_formatter(mpdates.DateFormatter( "%y%b" ) )
  fig.grid()
  if label :
    pyplot.ylabel( label )
  #fig.xaxis.set_major_formatter(mpdates.AutoDateFormatter())
  #pyplot.draw()

def RMplot( infile ) :
  t = convert( getcsvdata(infile,0))
  frac = numpy.array( getcsvdata(infile,5),dtype=float)
  pa = numpy.array( getcsvdata(infile,6),dtype=float)
  paerr = numpy.array( getcsvdata(infile,7),dtype=float)
  rm = numpy.array( getcsvdata(infile,8),dtype=float)
  rmerr = numpy.array( getcsvdata(infile,9),dtype=float)
  pyplot.clf()
  pyplot.ion()
  fig = pyplot.subplot(3,1,1)
  plot3( fig, t, rm, rmerr, 0, 15, label="RM (x 1.e5)")
  fig = pyplot.subplot(3,1,2)
  plot3( fig, t, pa, paerr, -90, 90, label="PA")
  fracerr = 0.001 * numpy.ones( len(frac) )
  fig = pyplot.subplot(3,1,3)
  plot3( fig, t, frac, fracerr, 0, .02, label="frac")
  pyplot.draw()
 
# convert from dateTime to Julian Day
# 2011 Jan 1 00:00:00 UT = jd 2455562.50
def dtToJD( dt ) :
  dtref = datetime.datetime.strptime("2011-jan-01-00:00","%Y-%b-%d-%H:%M")
  jdref = 2455562.50
  delt = dt - dtref
  return (jdref + delt.days + delt.seconds/86400.)
