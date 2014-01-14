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

