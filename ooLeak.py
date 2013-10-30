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
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages


class Leak:
  'leak object with data, legend, color'
 
  def __init__(self, file, antenna, legend, color ) :
    """read leakages from disk file, sort by frequency"""
    self.file = file
    self.ant = antenna
    self.legend = legend
    self.color = color
    solin = []
    self.f1 = []
    self.f2 = []
    self.DR = []
    self.DL = []
    self.avgchan = 0
    try :
      fin = open( self.file, "r" )
    except :
      print "... can't open file %s" % self.file
    else :
      print "... reading data from file %s" % self.file
      lastf2 = 0.
      for line in fin :
        a = line.split()
        if (line.startswith("# legend")) and (len(self.legend)==0) :
          self.legend = a[3]
        if line.startswith("# avgchan") :
          self.avgchan = int(a[3])
          #print "avgchan = %d" % self.avgchan
        if (len(a) > 9) and (line.startswith("C")) :		# new style Lk table, includes all antennas
          ant = int( a[0].strip("C") )
          if ant == self.ant :
            f1 = min(float(a[1]),float(a[2]))
            f2 = max(float(a[1]),float(a[2]))
            DR = float(a[3]) + 1j * float(a[4])
            DL = float(a[5]) + 1j * float(a[6])
            solin.append( [f1,f2,DR,DL] )
        if (not line.startswith("#")) and (len(a) > 10) :       # old style lk table, one antenna only
            f1 = min(float(a[0]),float(a[1]))
            f2 = max(float(a[0]),float(a[1]))
            DR = float(a[6]) + 1j * float(a[7])
            DL = float(a[8]) + 1j * float(a[9])
            chanfacts = a[12].split(",")
            self.avgchan = int(chanfacts[3])
            solin.append( [f1,f2,DR,DL] )
      fin.close()
    # it would be smarter to store values as self.sol = sorted..., but I am temporarily breaking
    #   up everything into f1,f2,DR,DL for compatibility with existing routines
      for s in sorted( solin, key = lambda(x) : x[1] ):	# store solutions in frequency order
        self.f1.append( s[0] )
        self.f2.append( s[1] )    
        self.DR.append( s[2] )
        self.DL.append( s[3] )
      #print self.f1
      #print self.f2
      #print self.DR
      #print self.DL
   
  def list(self) :
    return [self.ant, self.file, self.legend]

  def fminmax(self) :
    #fmin = self.sol[0][0]
    #fmax = self.sol[-1][1]
    fmin = self.f1[0]
    fmax = self.f2[-1]
    return [fmin,fmax] 
    

#  def dumpDR(self, fstart, fstop) :
#    if len(self.DR) <= 0 :
#      print "empty list"
#    else :
#      for f1,f2,DR in zip( self.f1, self.f2, self.DR ) :
#        if (f1 >= fstart) and (f1 <= fstop) and (f2 >=fstart) and (f2 <= fstop) :
#          print DR

#  def cDR(self, col) :
#    pylab.ion()
#    pylab.axis( [-.2, .2, -.2, .2] )
#    pylab.grid(True)
#    for DR in self.DR :
#      x = numpy.real(DR)
#      y = numpy.imag(DR)
#      pylab.plot(x,y,'bo',color=col,linestyle='solid',linewidth=1)
#    pylab.draw()

#  def cDL(self, col) :
#    x = []
#    y = []
#    pylab.ion()
#    pylab.axis( [-.2, .2, -.2, .2] )
#    pylab.grid(True)
#    for DL in self.DL :
#      x.append(numpy.real(DL))
#      y.append(numpy.imag(DL))
#    pylab.plot(x,y,color=col,linestyle='solid',linewidth=1)
#    pylab.draw()

#  def ampDR(self, col) :
#    pylab.ion()
#    pylab.axis( [self.fstart, self.fstop, 0.,.2] )
#    for f1,f2,DR in zip( self.f1, self.f2, self.DR ) :
#      amp = numpy.abs(DR)
#      #pylab.plot(f1,amp,'bo')
#      pylab.plot([f1,f2],[amp,mp],color=self.color,linestyle='solid',linewidth=1)
#    pylab.draw() 

# make bar graph of re and im vs freq
#  def riDL(self) :
#    pylab.ion()
#    pylab.axis( [self.fstart, self.fstop, -.1, .1] )
#    pylab.grid(True)
#    f2prev = 0.
#    for f1,f2,DL in zip( self.f1, self.f2, self.DL) :
#      y1 = numpy.real(DL)
#      y2 = numpy.imag(DL)
#      if f1 == f2prev :
#        pylab.plot( [f1,f1], [y1,y1prev] ,color='red',linestyle='solid',linewidth=1)
#        pylab.plot( [f1,f1], [y2,y2prev] ,color='blue',linestyle='solid',linewidth=1)
#      f2prev = f2
#      y1prev = y1
#      y2prev = y2
#      pylab.plot( [f1,f2], [y1,y1], color='red',linestyle='solid',linewidth=1)
#      pylab.plot( [f1,f2], [y2,y2], color='blue',linestyle='solid',linewidth=1)
#    pylab.draw()

  def plotComplex( self, p, fstart, fstop ) :
    """plot DR and DL on the complex plane"""
    first = True
    for f1,f2,DR,DL in zip( self.f1, self.f2, self.DR, self.DL ) :
      favg = 0.5 * (f1 + f2)
      if (favg > fstart) and (favg < fstop) :
        if first :     # plot dot and legend
          p.plot( DR.real, DR.imag, marker='o', color="red", markersize=5, \
             linestyle='solid', linewidth=1, label=self.legend )
          p.plot( DL.real, DL.imag, marker='D', color="blue", markersize=5 )
          first = False
        else :         
          if f1 == f2prev :			# connect with line
            p.plot( [xRprev,DR.real], [yRprev,DR.imag], marker='o', color="red", \
             markersize=5, linestyle='solid', linewidth=1 )
            p.plot( [xLprev,DL.real], [yLprev,DL.imag], marker='D', color="blue", \
             markersize=5, linestyle='dashed', linewidth=1 )
          else :                    # don't connect with line  
            p.plot( DR.real, DR.imag, marker='o', color="red", markersize=5 )
            p.plot( DL.real, DL.imag, marker='D', color="blue", markersize=5 )
        xRprev = DR.real
        yRprev = DR.imag
        xLprev = DL.real
        yLprev = DL.imag
        f2prev = f2
    p.legend( loc=0, prop={'size':10} )
         

  def panel(self, p, type, lk, fstart, fstop ) :
    """add to one panel of a plot; p = plot handle; type = amp,phs,complex; lk = DR or DL"""

    if lk == "DL" :  
      yc = self.DL
    else :
      yc = self.DR
    first = True
    f2prev = 0.

  # Go through the array freq by freq
    for f1,f2,ycomplex in zip( self.f1, self.f2, yc) :

        if (type == 'phs' ) : 
          y = numpy.angle(ycomplex, deg=True)
        else :
          y = numpy.abs(ycomplex)

      # Plot like histogram if these are > 240 MHz chunks
        if (abs(f1-f2) > .24) or (self.avgchan == 0) :
          if first :
            p.plot( [f1,f2], [y,y] , color=self.color, \
              linestyle='solid', linewidth=2, label=self.legend )
            first = False
          else :
            p.plot( [f1,f2], [y,y] , color=self.color, \
              linestyle='solid', linewidth=2 )
     
        else :  
          f = (f1+f2)/2.              # mean freq
        # Plot dots for phase
          if type == "phs" :
            p.plot( f, y, marker='o', color=self.color, markersize=3 )    # make dot
        # draw lines for amp
          else :
            if ((f1 > fstart) and (f1 < fstop)) or ((f2 > fstart) and (f2 < fstop)) :
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

  def loadAll(self, masterListFile ) :
    """Read in multiple leakage files to Plot object"""
    color = [ "black", "red", "blue", "green", "cyan", "magenta", "yellow" ]
    fin = open( masterListFile, "r" )
    ncolor = 0
    for line in fin :
      if len(line) > 0 :                    # skip blank lines
        truncate = line.find("#")           # skip anything after "#"
        a = line[0:truncate].split()
        if len(a) > 0 :
          for nant in range(1,16) :
            if a[0][-1] == "*" :            # old style lk files, one per antenna
              filename = a[0][0:-1] + str(nant)
            else :                          # new style Lk file, all antennas in 1 file
              filename = a[0]
            legend = ""
            if len(a) > 1 : legend = a[1]
            Plot.addLeak(self, filename, nant, legend, color[ncolor] )
              # addLeak opens disk file, reads data into Leak object
          ncolor = ncolor + 1
          if (ncolor > (len(color)-1) ) : ncolor = 0

  def list(self, antenna ) :
    n = 0
    for Leak in self.LeakList :
      [ant, file, legend] = Leak.list()
      if ant == antenna :
        n = n + 1
        print n, self.plotList[n-1], file, legend

#  def amp(self, antenna, f1=0., f2=0., amax=.25) :
#    [f1, f2] = Plot.xlimits( self )
#    pylab.clf()
#    pylab.ion()
#    for Leak in self.LeakList :
#      Leak.ap( 'amp', antenna, f1, f2, 0., amax )

#  def phs(self, antenna, f1=0., f2=0.) :
#    [f1, f2] = Plot.xlimits( self )
#    pylab.clf()
#    pylab.ion()
#    for Leak in self.LeakList :
#      Leak.ap( 'phs', antenna, f1, f2, -180, 180. )

  def xlimits(self) :
    f1 = 300.
    f2 = 0.
    for Leak in self.LeakList :
      [fmin,fmax] = Leak.fminmax()
      if fmin < f1 : f1 = fmin
      if fmax > f2 : f2 = fmax
    f1 = f1 - 0.05*(f2-f1) 
    f2 = f2 + 0.05*(f2-f1) 
    return [f1, f2] 

  def replot(self, fstart, fstop ) :
    pylab.ion()
    pylab.axis( [fstart, fstop, 0, .25] )
    pylab.grid(True)
   
  def clear(self) :
    pylab.clf()

  def ampAll(self, f1=0., f2=0., type="amp", amax=.25) :
    pyplot.ioff()
    pp = PdfPages( 'AmpLeaks.pdf' )
    if f1 == 0. :
      [f1, f2] = Plot.xlimits( self )          # default is to find freq limits in the data
      print f1,f2
    ymin = 0.
    ymax = amax
    if type == "phs" :
      ymin = -180.
      ymax = 180.
    for ant in range(1,16) :
      pyplot.clf()
      pL = pyplot.subplot(2,1,1)    # DL in upper panel
      pL.axis( [f1, f2, ymin, ymax] )
      pL.grid(True)
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DL for antenna %d" % ant
          Leak.panel(pL, type, "DL", f1, f2 )
          pyplot.title("C%d DL" % ant)
      pR = pyplot.subplot(2,1,2)    # DR in lower panel
      pR.axis( [f1, f2, ymin, ymax] )
      pR.grid(True)
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DR for antenna %d" % ant
          Leak.panel(pR, type, "DR", f1, f2 )
          pyplot.title("C%d DR" % ant)
      pyplot.savefig( pp, format='pdf' )
    pp.close()

# for complex plot, f1 and f2 could in principle be used to select range of points
#    but this is not implemented yet

  def complexAll(self, f1=0., f2=0., type="complex", amax=.16) :
    pyplot.ioff()
    pp = PdfPages( 'ComplexLeaks.pdf' )
    if f1 == 0. :
      [f1, f2] = Plot.xlimits( self )          # default is to find freq limits in the data
      print f1,f2
    ymin = -1.*amax
    ymax = amax
    for ant in range(1,16) :
      pyplot.clf()
      p = pyplot.subplot(1,1,1, aspect='equal')    # DL,DR in one panel
      p.axis( [ymin, ymax, ymin, ymax] )
      p.grid(True)
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DR and DL for antenna %d" % ant
          Leak.plotComplex( p, f1, f2 ) 
          #Leak.panel(p, type, "DL", f1, f2 )
          #Leak.panel(p, type, "DR", f1, f2 )
      pyplot.title("C%d DR (circles, solid) and DL (diamonds, dashed)" % ant)
      pyplot.savefig( pp, format='pdf' )
    pp.close()
