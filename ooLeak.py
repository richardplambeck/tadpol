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
 
  def __init__(self, file, antenna, legend, color, marker ) :
    """read leakages from disk file, sort by frequency"""
    self.file = file
    self.ant = antenna
    self.legend = legend
    self.color = color
    self.marker = marker
    solin = []
    self.f1 = []
    self.f2 = []
    self.DR = []
    self.DL = []
    self.avgchan = "0"
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
          self.avgchan = a[3]
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
            self.avgchan = chanfacts[3]
            solin.append( [f1,f2,DR,DL] )
      fin.close()
    # it would be smarter to store values as self.sol = sorted..., but I am temporarily breaking
    #   up everything into f1,f2,DR,DL for compatibility with existing routines
      for s in sorted( solin, key = lambda(x) : x[1] ):	# store solutions in frequency order
        self.f1.append( s[0] )
        self.f2.append( s[1] )    
        self.DR.append( s[2] )
        self.DL.append( s[3] )
   
  def list(self) :
    return [self.ant, self.file, self.legend]

  def fminmax(self) :
    fmin = self.f1[0]
    fmax = self.f2[-1]
    return [fmin,fmax] 
    
  def plotComplex( self, p, fstart, fstop ) :
    """plot DR and DL on the complex plane"""
    first = True
    f2prev = 0.
    DRmark = "."
    DLmark = "D"
    for f1,f2,DR,DL in zip( self.f1, self.f2, self.DR, self.DL ) :
      msize = 5
      #if (f2-f1) > 0.24 :
      #  msize = 8
      favg = 0.5 * (f1 + f2)
      if (favg > fstart) and (favg < fstop) and (numpy.abs(DR) != 0.) :
        if first :     # plot dot and legend
          p.plot( DR.real, DR.imag, marker=DRmark, color=self.color, markersize=msize, \
             label=self.legend )
          p.plot( DL.real, DL.imag, marker=DLmark, color=self.color, markersize=msize )
          first = False
        else :         
          #if f1 == f2prev :			# connect with line
          if True :
            p.plot( [xRprev,DR.real], [yRprev,DR.imag], marker=DRmark, color=self.color, \
             markersize=msize, linestyle='solid', linewidth=1 )
            p.plot( [xLprev,DL.real], [yLprev,DL.imag], marker=DLmark, color=self.color, \
             markersize=msize, linestyle='dashed', linewidth=1 )
          else :                    # don't connect with line  
            p.plot( DR.real, DR.imag, marker=DRmark, color=self.color, markersize=msize )
            p.plot( DL.real, DL.imag, marker=DLmark, color=self.color, markersize=msize )
        xRprev = DR.real
        yRprev = DR.imag
        xLprev = DL.real
        yLprev = DL.imag
        f2prev = f2
    p.legend( loc=0, prop={'size':6}, numpoints=1 )
         

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
        if (abs(f1-f2) > .24) or (self.avgchan == "0") :
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
  def addLeak(self, file, antenna, legend, color, marker ) :
    newLeak = Leak( file, antenna, legend, color, marker )
    self.LeakList.append( newLeak )
    self.plotList.append( True )

  def loadAll(self, masterListFile ) :
    """Read in multiple leakage files to Plot object"""
    color = [ "black", "red", "blue", "green", "cyan", "magenta", "yellow" ]
    marker = [ "o", "D", "v", "^", "s", "h", "d" ]
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
            Plot.addLeak(self, filename, nant, legend, color[ncolor], marker[ncolor] )
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

  def complexAll(self, f1=0., f2=0., amax=.16, nrows=1, ncols=1) :
    pyplot.ioff()
    pp = PdfPages( 'ComplexLeaks.pdf' )
    scale = 10./math.sqrt(ncols*nrows)
    if f1 == 0. :
      [f1, f2] = Plot.xlimits( self )          # default is to find freq limits in the data
    ymin = -1.*amax
    ymax = amax
    npanel = 0
    for ant in range(1,16) :
      npanel = npanel + 1
      if npanel > nrows * ncols :
        npanel = 1
        pyplot.clf()
      p = pyplot.subplot(nrows, ncols, npanel, aspect='equal')    # DL,DR in one panel
      p.tick_params( axis='both', which='major', labelsize=scale )
      p.axis( [ymin, ymax, ymin, ymax] )
      p.grid(True)
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DR and DL for antenna %d" % ant
          Leak.plotComplex( p, f1, f2 ) 
      pyplot.title("C%d DR (circles, solid) and DL (diamonds, dashed)" % ant, fontdict={'fontsize': scale})
      if (npanel == nrows*ncols) or (ant == 15) :
        pyplot.savefig( pp, format='pdf' )
    pp.close()

# first file in LeakList is Lktemplate
# for each freq interval in Lktemplate that is within range f1-f2, find all other entries in LeakList
#   which match within 10% of the freq range; append these to list; form avg of list; plot scatter
#   relative to the avg, 15 plots per page; compute rms

  def complexCmp(self, fmin=0., fmax=300., amax=.05, nrows=4, ncols=4) :
    Lktemplate = self.LeakList[0]
    DRmark = "o"
    DLmark = "x"
    pyplot.ioff()
    pp = PdfPages( 'ComplexCmp.pdf' )
    for f1ref,f2ref in zip ( Lktemplate.f1, Lktemplate.f2 ) :
      if (f1ref > fmin) and (f2ref < fmax) :
        finterval = f2ref-f1ref
        print f1ref,f2ref,finterval
        npanel = 0
        for ant in range(1,16) :
          DRlist = []
          DLlist = []
          leglist = []
          collist = []
          for Lk in self.LeakList :
            if Lk.ant == ant :
              for f1,f2,DR,DL in zip( Lk.f1, Lk.f2, Lk.DR, Lk.DL ) :
                if (abs( (f1-f1ref)/finterval ) < 0.1) and (abs( (f2-f2ref)/finterval) < 0.1) and (abs(DR) > 0.) and (abs(DL) > 0.) :
                  DRlist.append( DR )
                  DLlist.append( DL )
                  leglist.append(Lk.legend)
                  collist.append(Lk.color)
                  print ant,DR,DL,Lk.legend
          DRmean = numpy.mean(DRlist)
          DLmean = numpy.mean(DLlist)
          print ant,DRmean,DLmean
          npanel= npanel + 1
          p = pyplot.subplot(nrows, ncols, npanel, aspect='equal')    # DL,DR in one panel
          p.tick_params( axis='both', which='major', labelsize=5 )
          p.axis( [-1.*amax, amax, -1.*amax, amax] )
          p.grid(True)
          p.text(.93,.81,"C%d" % ant, horizontalalignment='right',transform=p.transAxes, size=7)
          for DR,legend,color in zip( DRlist,leglist,collist ) :
            p.plot( (DR-DRmean).real, (DR-DRmean).imag, marker=DRmark, color=color, markersize=5, \
              label=legend )
          for DL,color in zip( DLlist,collist) :
            p.plot( (DL-DLmean).real, (DL-DLmean).imag, marker=DLmark, color=color, markersize=5 )
          if ant == 13 :
            p.legend(bbox_to_anchor=(5.7,0.9), prop={'size':5})
            p.text(5.5,0.92,"(%.3f,%.3f)" % (f1ref,f2ref), horizontalalignment='center',transform=p.transAxes, size=7)
        pyplot.savefig( pp, format='pdf' )
        pyplot.clf()
    pp.close()
         
          
