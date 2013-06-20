# ooLeak.py
#
#import subarrayCommands as SAC
import math
import time
#import device
import cmath
import numpy
import pylab
import sys

# a 'Leak' object contains the leakage solution for 1 antenna on 1 night
# ideally, this solution has a particular freq resolution, although nothing
#   prevents several freq resolutions, or even duplicate frequency bands,
#   in one object
# typically the 'legend' holds the date and src

class Leak:
  'leak object with data, legend, color'
 
  def __init__(self, file, antenna, legend, color ) :
    self.file = file
    self.ant = antenna
    self.legend = legend
    self.color = color
    self.f1 = []
    self.f2 = []
    self.DR = []
    self.DL = [] 
    try :
      fin = open( self.file, "r" )
    except :
      print "... can't open file %s" % self.file
    else :
      print "... reading data from file %s" % self.file
      lastf2 = 0.
      for line in fin :
        a = line.split()
        if (not line.startswith("#")) and (len(a) > 10) :
          self.f1.append(float(a[0]))
          self.f2.append(float(a[1]))
          self.DR.append(float(a[6]) + 1j * float(a[7]))
          self.DL.append(float(a[8]) + 1j * float(a[9]))
      fin.close()

  def list(self) :
    return [self.ant, self.file, self.legend]

  def DR(self) :
    return self.DR

  def minmax(self) :
    fmin = 1000.
    fmax = 0.
    for f1,f2 in zip( self.f1, self.f2 ) :
      if f1 < fmin : fmin = f1
      if f2 < fmin : fmin = f2 
      if f1 > fmax : fmax = f1
      if f2 > fmax : fmax = f2
    return [fmin,fmax] 
    

  def dumpDR(self, fstart, fstop) :
    if len(self.DR) <= 0 :
      print "empty list"
    else :
      for f1,f2,DR in zip( self.f1, self.f2, self.DR ) :
        if (f1 >= fstart) and (f1 <= fstop) and (f2 >=fstart) and (f2 <= fstop) :
          print DR

  def cDR(self, col) :
    pylab.ion()
    pylab.axis( [-.2, .2, -.2, .2] )
    pylab.grid(True)
    for DR in self.DR :
      x = numpy.real(DR)
      y = numpy.imag(DR)
      pylab.plot(x,y,'bo',color=col,linestyle='solid',linewidth=1)
    pylab.draw()

  def cDL(self, col) :
    x = []
    y = []
    pylab.ion()
    pylab.axis( [-.2, .2, -.2, .2] )
    pylab.grid(True)
    for DL in self.DL :
      x.append(numpy.real(DL))
      y.append(numpy.imag(DL))
    pylab.plot(x,y,color=col,linestyle='solid',linewidth=1)
    pylab.draw()

  def ampDR(self, col) :
    pylab.ion()
    pylab.axis( [self.fstart, self.fstop, 0.,.2] )
    for f1,f2,DR in zip( self.f1, self.f2, self.DR ) :
      amp = numpy.abs(DR)
      #pylab.plot(f1,amp,'bo')
      pylab.plot([f1,f2],[amp,amp],color=self.color,linestyle='solid',linewidth=1)
    pylab.draw() 

  def setrange(self, fstart, fstop ) :
    self.fstart = fstart
    self.fstop = fstop

# make bar graph of re and im vs freq
  def riDL(self) :
    pylab.ion()
    pylab.axis( [self.fstart, self.fstop, -.1, .1] )
    pylab.grid(True)
    f2prev = 0.
    for f1,f2,DL in zip( self.f1, self.f2, self.DL) :
      y1 = numpy.real(DL)
      y2 = numpy.imag(DL)
      if f1 == f2prev :
        pylab.plot( [f1,f1], [y1,y1prev] ,color='red',linestyle='solid',linewidth=1)
        pylab.plot( [f1,f1], [y2,y2prev] ,color='blue',linestyle='solid',linewidth=1)
      f2prev = f2
      y1prev = y1
      y2prev = y2
      pylab.plot( [f1,f2], [y1,y1], color='red',linestyle='solid',linewidth=1)
      pylab.plot( [f1,f2], [y2,y2], color='blue',linestyle='solid',linewidth=1)
    pylab.draw()

  def panel(self, p, type, lk, fstart, fstop, ymin, ymax ) :
    """one plot panel; type = amp,phs"""
    p.axis( [fstart, fstop, ymin, ymax] )
    p.grid(True)
    if lk == "DL" :  
      yc = self.DL
    else :
      yc = self.DR
    first = True
    f2prev = 0.
    for f1,f2,ycomplex in zip( self.f1, self.f2, yc) :
      if (type == 'phs' ) : 
        y = numpy.angle(ycomplex, deg=True)
      else :
        y = numpy.abs(ycomplex)
      if abs(f1-f2) > .24 :
        if first :
          p.plot( [f1,f2], [y,y] , color=self.color, \
            linestyle='solid', linewidth=2, label=self.legend )
          first = False
        else :
          p.plot( [f1,f2], [y,y] , color=self.color, \
            linestyle='solid', linewidth=2 )
      else :  
        f = (f1+f2)/2.              # mean freq
        if ((f1 > fstart) and (f1 < fstop)) or ((f2 > fstart) and (f2 < fstop)) :
          #pL.plot( f, y, marker='.', color=self.color )    # make dot
          if f1 == f2prev :
            if first :
              p.plot( [fprev,f], [yprev,y] , color=self.color, \
                linestyle='solid', linewidth=1, label=self.legend )
              first = False
            else :
              p.plot( [fprev,f], [yprev,y] , color=self.color, \
                linestyle='solid', linewidth=1 )
          fprev = f
          f2prev = f2
          yprev = y
    p.legend( loc=0, prop={'size':10} )
    
  def ap(self, type, antenna, fstart, fstop, ymin, ymax )  :
    """plot amplitude or phase of DL and DR for one antenna"""
    if self.ant != antenna :
      return 
    pL = pylab.subplot(2,1,1)    # DL in upper panel
    Leak.panel(self, pL, type, "DL", fstart, fstop, ymin, ymax )
    pylab.title("C%d DL" % antenna)
    pR = pylab.subplot(2,1,2)    # DR in lower panel
    Leak.panel(self, pR, type, "DR", fstart, fstop, ymin, ymax )
    pylab.title("C%d DR" % antenna)
    print "begin draw"
    pylab.draw()
    

# ====================================================================================#

class Plot:
  """Plot object contains one or more Leak objects"""
 
  # --- initialize plot object; optionally, load leakage data from list
  def __init__(self, masterList ) :
    self.LeakList = []
    self.plotList = []
    self.fstart = 200.
    self.fstop = 270.
    if masterList == None :
      print "creating blank plot object"
    else :
      print "loading leakage objects from file %s" % masterList
      Plot.loadAll(self, masterListFile=masterList )

  # --- add one Leak to Plot object ---
  def addLeak(self, file, antenna, legend, color ) :
    newLeak = Leak( file, antenna, legend, color )
    self.LeakList.append( newLeak )
    self.plotList.append( True )

  # --- load multiple Leaks to Plot object ---
  def loadAll(self, masterListFile ) :
    color = [ "black", "red", "blue", "green", "cyan", "magenta", "yellow" ]
    fin = open( masterListFile, "r" )
    ncolor = 0
    for line in fin :
      if not line.startswith("#") :
        a = line.split()
        if len(a) > 0 :
          for nant in range(1,16) :
            filename = a[0] + str(nant)
            Plot.addLeak(self, filename, nant, a[1], color[ncolor] )
          ncolor = ncolor + 1
          if (ncolor > (len(color)-1) ) : ncolor = 0

  def list(self, antenna ) :
    n = 0
    for Leak in self.LeakList :
      [ant, file, legend] = Leak.list()
      if ant == antenna :
        n = n + 1
        print n, self.plotList[n-1], file, legend

  def amp(self, antenna, f1=0., f2=0., amax=.25) :
    [f1, f2] = Plot.xlimits( self, antenna, f1, f2 )
    pylab.clf()
    pylab.ion()
    for Leak in self.LeakList :
      Leak.ap( 'amp', antenna, f1, f2, 0., amax )

  def phs(self, antenna, f1=0., f2=0.) :
    [f1, f2] = Plot.xlimits( self, antenna, f1, f2 )
    pylab.clf()
    pylab.ion()
    for Leak in self.LeakList :
      Leak.ap( 'phs', antenna, f1, f2, -180, 180. )

  def xlimits(self, antenna, f1, f2) :
    if (f1 == 0.) :
      f1 = 300.
      for Leak in self.LeakList :
        if Leak.ant == antenna :
          [fmin,fmax] = Leak.minmax()
          if fmin < f1 : f1 = fmin
          if fmax > f2 : f2 = fmax
      #print "%.3f, %.3f" % (fstart,fstop)
      f1 = f1 - 0.05*(f2-f1) 
      f2 = f2 + 0.05*(f2-f1) 
      #print "%.3f, %.3f" % (fstart,fstop)
    return [f1, f2] 

  def replot(self, fstart, fstop ) :
    pylab.ion()
    pylab.axis( [fstart, fstop, 0, .25] )
    pylab.grid(True)
   
  def clear(self) :
    pylab.clf()

  def DRampAll(self, f1=0., f2=0., amax=.25) :
    pylab.clf()
    pylab.ion()
    [f1, f2] = Plot.xlimits( self, 1, f1, f2 )
    for ant in range(1,16) :
      pAnt = pylab.subplot(5,3,ant)   
      pylab.title("C%d DR" % ant)
      for Leak in self.LeakList :
        if Leak.ant == ant :
          Leak.panel( pAnt, "amp", "DR", f1, f2, 0., amax )
          pylab.draw()