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
import vCal
import string
import leakSolve
import paPlot
import pickle
import datetime
import scipy.optimize
import dateutil
import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
import matplotlib.dates as mpdates
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import LogFormatterExponent

     
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


def getcsvdata(infile,col,type=None) :
  s = []
  fin = open(infile,"r")
  for line in fin:
    if not line.startswith("#"):
      a = line.split(',')
      if type == "Float" :
        s.append( float(a[col]) )
      else :
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
def dt2JD( dt ) :
  dtref = datetime.datetime.strptime("2011-jan-01-00:00","%Y-%b-%d-%H:%M")
  jdref = 2455562.50
  delt = dt - dtref
  return (jdref + delt.days + delt.seconds/86400.)

# find distance between Jupiter and Ganymede, say
def separation( ref, source, outfile ) :
  fout = open( outfile, "w" )
  for mjdiff in numpy.arange(0.,40.,.1) :
    p= subprocess.Popen( ( shlex.split('distance ref=%s source=%s mjd=+%.2f'  \
       % (ref, source, deltamjd) ) ), \
       stdout=subprocess.PIPE,stdin=subprocess.PIPE,stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    pout.write("%s\n" % result)
  fout.close()

def plot4( t, s) :
  pyplot.clf()
  pyplot.ion()
  fig = pyplot.subplot(1,1,1)
  tt = mpdates.date2num(t)
  pyplot.ylim(0,10.)
  fig.plot_date( tt, s, "r." )
  #fig.xaxis.set_major_locator(mpdates.MinuteLocator(interval=10))
  fig.grid()
  pyplot.ylabel("angular offset (arcmin)")
  pyplot.xlabel("date")
  #fig.xaxis.set_major_locator(mpdates.AutoDateLocator())
  #fig.xaxis.set_major_formatter(mpdates.AutoDateFormatter())
  fig.xaxis.set_major_formatter(mpdates.DateFormatter( "%m/%d" ) )
  pyplot.draw()
  pyplot.show()

def plotsep( infile ) :
  t = []
  s = []
  fin = open( infile, "r" )
  for line in fin:
    if line.startswith("2015") :
      a = line.split()
      elev = float(a[2])
      if elev > 10. :
        t.append(datetime.datetime.strptime( a[0][0:16], "%Y-%m-%dT%H:%M") )
        s.append(60.*float(a[5]))
  plot4( t,s )
  fin.close()

def plotHA( infile="List2" ) :
  times = numpy.array( [0,0], dtype=float )
  pyplot.clf()
  pyplot.ion()
  fig = pyplot.subplot(1,1,1)
  pyplot.ylim( 0., 25. )
  pyplot.xlim( 0., 24. )
  fig.grid()
  yoff = [1.,1.]
  fin = open( infile, "r" )
  lastSrc = None
  for line in fin :
    if string.find( line, "NEVER" ) < 0 :
       a = line.split()
       color = 'red'
       if float(a[9]) > 15. : color = 'green'
       if float(a[9]) > 20. : color = 'blue'
       if lastSrc != a[0] :
         yoff[0] = yoff[0] + 1.
         yoff[1] = yoff[1] + 1.
         fig.text( vCal.dechrs(a[5])-.2, yoff[0], a[0], verticalalignment="center", \
             horizontalalignment="right", fontsize=10 )
         lastSrc = a[0]
       times[0] = vCal.dechrs(a[5])
       times[1] = vCal.dechrs(a[6])
       
       if times[1] < times[0] :
         times[0] = 0.
         fig.plot( times, yoff, linewidth=5, color=color )
         times[0] = vCal.dechrs(a[5])
         times[1] = 24.
       fig.plot( times, yoff, linewidth=5, color=color )
  fin.close()
  pyplot.draw()

# used by vCal to plot smoothed gains
# g1 is original gain
# g2 is new gain
def plotGains( t, g1, g2, plotName ) :
  nrows = 2
  ncols = 2
  pyplot.clf()
  pyplot.ioff()
  pp = PdfPages( plotName )
  nplot = 1
  for nant in range(1,16) :
    if nplot > nrows*ncols :
      pyplot.suptitle( plotName, size=8 )
      pyplot.savefig( pp, format='pdf' )
      nplot = 1
      pyplot.clf()
    fig=pyplot.subplot( nrows, ncols, nplot )
    fig.tick_params( axis='both', which='major', labelsize=8)
    if nplot > 2 :
      pyplot.xlabel( 'UT hrs', fontsize=8 )
    if nplot == 1 or nplot == 3 :
      pyplot.ylabel( 'gain correction (sqrt Jy/K)', fontsize=8 )
    pyplot.ylim( -0.1, 3.1)
    fig.plot( t, g1[:,[nant-1]], "r,")
    fig.plot( t, g2[:,[nant-1]], "b," )
    # fig.plot( t, g3[ nant-1 ], "-" )
    fig.text( .05, .9, "C%d" % nant, transform=fig.transAxes )
    nplot = nplot + 1
  #pyplot.draw()
  if (nplot > 1) :
    pyplot.suptitle( plotName, size=8 )
    pyplot.savefig( pp, format='pdf' )
  pp.close()

def plotSEFD( infile, label ) :
    fin = open( infile, "r" ) 
    pyplot.clf()
    pyplot.ion()
    fig = pyplot.subplot(1,1,1)
    for line in fin :
      if not line.startswith('#') :
        a = line.split()
        color = 'black'
        if (a[2] == "SGRA" ) : color = 'red'
        if (a[2] == "M87" ) :  color = 'blue'
        t = datetime.datetime.strptime( a[3], "%H:%M:%S") 
        tt = mpdates.date2num(t)
        Cr = float(a[7]) 
        Cp = float(a[8]) 
        Cplow = float(a[9]) 
        Cphigh = float(a[10]) 
        fig.plot_date( tt, Cp, marker="o", color=color, linewidth=3 )
        fig.plot_date( tt, Cr, marker="s", fillstyle='full', alpha=0.5, color=color )
        fig.plot_date( [tt,tt], [Cplow,Cphigh], "-", linewidth=2, color=color )
    fin.close()
    fig.xaxis.set_major_locator(mpdates.HourLocator())
    fig.xaxis.set_major_formatter(mpdates.DateFormatter( "%H" ) )
    fig.set_yscale("log")
    tmin = mpdates.date2num( datetime.datetime.strptime( "00:00", "%H:%M") )
    tmax = mpdates.date2num( datetime.datetime.strptime( "19:00", "%H:%M") )
    fig.set_xlim( tmin-.02,tmax )
    fig.set_ylim( 1.e3,1.e6 )
    fig.text( .5, .94, label, verticalalignment="center", transform=fig.transAxes, \
      horizontalalignment="center", fontsize=18 )
    
    fig.grid()
    pyplot.xlabel( "UT (hrs)" )
    pyplot.ylabel( "SEFD (Jy)" )
    pyplot.show()
   
# for PB2 optics tube
keyList = [ [2,"50K coldhead","slateblue","-"],
            [5,"50K link cold end","sienna","-"],
            [7,"50K link hot end","red","-"],
            [8,"50K link middle","gold","-"],
           [1,"4K coldhead","lime","-"],
           [3,"4K link","springgreen","--"],]

# routines for BlackDewar 
#keyList = [ [1,"interhead","lime","-"],
#            [2,"ultra 2 U04577","slateblue","--"],
#            [3,"ultra 3 31275","springgreen","--"],
#            [4,"ultrahead","blue","-"],
#            [7,"4He exchanger","red","-"],
#            [8,"4He film burner","gold","-"],
#            [9,"4K shield sensor 4","magenta","--"],
#            [13,"PT 4K head D07697","sienna","-"] ]

#            [14,"PT 50K head","darkred","-"]]
#            [5,"4He pump","red","--"],
#            [6,"4He switch","red",":"],
#            [10,"3He IC switch","lime",":"],
#            [11,"3He UC pump","blue","--"],
#            [12,"3He UC switch","blue",":"],
#            [15,"sensor 3","fuchsia","-"], 
#            [16,"sensor 2","indigo","-"] ]

#keyList = [ [1,"interhead","lime","-"],
#            [2,"TowerTop","slateblue","--"],
#            [3,"TowerMid","springgreen","--"],
#            [4,"ultrahead","blue","-"],
#            [7,"4He exchanger","red","-"],
#            [14,"TowerBase","magenta","-"] ]

#keyList = [ [1,"interhead","lime","-"],
#            [2,"ultrahead","blue","-"],
#            [3,"4He pump","red","--"],
#            [4,"4He pump switch","red",":"],
#            [5,"4He exchanger","red","-"],
#            [6,"4He film burner","tomato","-"],
#            [8,"3He interpump switch","lime",":"],
#            [9,"3He ultrapump","blue","--"],
#            [10,"3He ultrapump switch","blue",":"],
#            [11,"PT 4K head","peru","-"],
#            [12,"tower 2K","maroon","-"] ]
#            [7,"3He interpump","lime","--"],

# reads cryoLog text file, returns array of datetimes and masked array of temperatures
# only columns listed in keyList will be returned
# BD=True for BlackDewar format, False for PB2OT format
def getBDdata( infile, keyList=keyList, BD=False ) :
  t = []
  s = []
  First = True
  t0 = 0.
  fin = open(infile,"r")
  for line in fin:
    if BD :
      if not (line.startswith("#") or line.startswith("Time") ) :
        a = line.split()
        t.append(datetime.datetime.strptime(a[0],"%Y%m%d%H%M%S"))
        sList = []
        for col in range(1,len(a)) :
          sList.append(float(a[col]))
        #print t,sList
        s.append(sList)
    else :
      if not (line.startswith("Date") or line.startswith("_") or line.startswith("Time") ) :
        a = line.split(',')
        if First :
          t0 = float(a[0])
          First = False
        t.append( float(a[0]) - t0 )
        sList = []
        for col in range(1,len(a)) :
          sList.append(float(a[col]))
        #print t[-1],sList
        s.append(sList)
  fin.close()
  return t,numpy.ma.masked_less_equal( numpy.array(s), 0. )

def copyValues( t, s, ncolumn, smin, smax ) :
    tplot = []
    splot = []
    for i in range(0,len(t)) :
      if (s[i,ncolumn-1] >= smin) and (s[i,ncolumn-1] <= smax) :
        tplot.append( t[i] )
        splot.append( s[i,ncolumn-1] )
    return tplot,splot

def plotBD( infile, keyList=keyList, plotTitle=None, time1=0., time2=0., temperature1=0.2, temperature2=300., BD=False ) :
    pyplot.ioff()
    pp = PdfPages( "cryoLog.pdf" )
    fig, ax = pyplot.subplots( figsize=(11,8) )
    t,s = getBDdata( infile, keyList=keyList )
    print s[:,0]
    print s[:,1]
    print s.shape
    nc = 0 
    for [column,name,color,linetype] in keyList :
      print column, name, color
      tplot,splot = copyValues( t, s, column, 0., 300. )
      if BD :
        ax.plot_date( tplot, splot, fmt='-', color=color, linewidth=2, linestyle=linetype, label=name )
      else :
        ax.plot( tplot, splot, '-', color=color, linewidth=2, linestyle=linetype, label=name )
      nc = nc + 1
    if BD :
      ax.xaxis.set_major_locator(mpdates.HourLocator())
      ax.xaxis.set_major_formatter(mpdates.DateFormatter( "%H" ) )
      datelabel = t[0].strftime("%Y %B %d %H:%M")
      if time1 != 0. or time2 != 0. : 
        time1 = datetime.datetime.strptime(time1,"%Y%m%d%H%M%S")
        datelabel = time1.strftime("%Y %B %d %H:%M")
        time2 = datetime.datetime.strptime(time2,"%Y%m%d%H%M%S")
        ax.set_xlim( time1,time2 )
      pyplot.xlabel( "PDT beginning %s" % datelabel )
    else :
      pyplot.xlabel( "elapsed time (secs)" )
    ax.tick_params( axis='y', which='minor' )
    #ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    #ax.set_yscale("log")
    #ax.yaxis.set_minor_formatter( LogFormatterExponent )
    if plotTitle:
      fig.text( .5, .92, plotTitle, verticalalignment="center", transform=ax.transAxes, \
        horizontalalignment="center", fontsize=14 )
    
    #ax.set_ylim( temperature1,temperature2 )
    ax.grid( which='both' )
    pyplot.ylabel( "temperature (K)" )
    pyplot.legend( loc="best", prop={'size':10} )
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()
   
# derive calibration curve for a sensor
def calCurve( infile1, infile2, Tcol, Rcol ) :
    fin1 = open( infile1, "r" )
    fin2 = open( infile2, "r" )
    fout = open( "junk", "w" )
    nr1 = 5 + 10*Tcol 
    ns1 = 5 + 10*Rcol
    nl = 0
    for line1,line2 in zip(fin1,fin2) :
      if not line1.startswith("#") :
        nl = nl + 1
        print line1[0:14], line2[0:14], line1[nr1:nr1+10], line2[ns1:ns1+10] 
        fout.write("%5d %s %s %s %s\n" % (nl,line1[0:14],line2[0:14],line1[nr1:nr1+10],line2[ns1:ns1+10]) )

# for PB2 optics tube
keyList1 = [ [2,"50K coldhead","slateblue","-"],
            [5,"50K link cold end","sienna","-"],
            [7,"50K link hot end","red","-"],
            [8,"50K link middle","gold","-"]]
keyList2 = [ [1,"4K coldhead","blue","-"],
           [3,"4K link","lime","-"],]
    
def plotOT( infile, keyList1=keyList1, keyList2=keyList2, plotTitle=None, time1=0., time2=0., temperature1=0.2, temperature2=300., BD=False ) :
    pyplot.ioff()
    pp = PdfPages( "cryoLog.pdf" )
    fig, (ax1,ax2) = pyplot.subplots( 2, 1, figsize=(11,8), sharex='all', squeeze=True )
    t,s = getBDdata( infile, keyList=keyList1 )
    print s[:,0]
    print s[:,1]
    print s.shape
    nc = 0 
    for [column,name,color,linetype] in keyList1 :
      print column, name, color
      tplot,splot = copyValues( t, s, column, 0., 300. )
      ax1.plot( tplot, splot, '-', color=color, linewidth=2, linestyle=linetype, label=name )
      nc = nc + 1
    ax1.tick_params( axis='y', which='minor' )
    if plotTitle:
      fig.text( .5, .92, plotTitle, verticalalignment="center", transform=ax.transAxes, \
        horizontalalignment="center", fontsize=14 )
    ax1.grid( which='both' )
    ax1.set_ylabel( "temperature (K)" )
    ax1.legend( loc="best", prop={'size':10} )

    nc = 0 
    for [column,name,color,linetype] in keyList2 :
      print column, name, color
      tplot,splot = copyValues( t, s, column, 0., 300. )
      ax2.plot( tplot, splot, '-', color=color, linewidth=1, linestyle=linetype, label=name )
      nc = nc + 1
      pyplot.xlabel( "elapsed time (secs)" )
    ax2.tick_params( axis='y', which='minor' )
    ax2.grid( which='both' )
    pyplot.xlabel( "elapsed time (secs)" )

    pyplot.ylabel( "temperature (K)" )
    pyplot.legend( loc="best", prop={'size':10} )
    pyplot.savefig( pp, format='pdf' )
    pp.close()
    pyplot.show()
   
