# ooLeak.py

#
# this is a collection of routines to plot, list, compare, or average leakages
# unfortunately I chose to make it object-oriented, which turned out not to be particularly useful
#
# there are 2 objects:
# - a Leak object contains the frequency-dependent leakages for a single antenna, derived from
#     observations of a particular source, with a particular channel averaging interval
# - a LkSet object is a collection of multiple Leak objects
#
# the actual leakage solutions are computed by leakSolve, and are written to disk files with 
#   names like, e.g., Lka.21mar2013; an older format with one antenna per file, e.g., lk13, also
#   may be read in
#
# to get away from the object-oriented stuff, I have some wrapper routines that do what I usually
#   want to do: plotAllAmps, etc; these are at the beginning 


import math
import time
import cmath
import numpy
import pylab
import sys
import RM
import pickle
import random
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle

allAnts = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]   # default antList

# ----------------------------------------------------------------------------------------------------- #
# the Leak object - leakages vs frequency for one antenna
# when given a standard Lk file, opens it, strips out only the data for specfied antenna

class Leak:
  'leak object with data, legend, color'
 
  def __init__(self, file, antenna, legend, color, marker ) :
    """read leakages from disk file, sort by frequency"""
    self.file = file        # name of file in quotes
    self.ant = antenna      # integer
    self.legend = legend
    self.color = color
    self.marker = marker
    self.selectStr = ""
    solin = []	             # just a sequential list
    self.f1 = []
    self.f2 = []
    self.DR = []
    self.DL = []
    self.Qpercent = []
    self.Upercent = []
    self.lineStr = []
    self.avgchan = "0"
    try :
      fin = open( self.file, "r" )
    except :
      print "... can't open file %s" % self.file
    else :
      print "... reading leakages for antenna %d from file %s" % (antenna,self.file)
      lastf2 = 0.
      for line in fin :
        a = line.split()
        if (line.startswith("# legend")) and (len(self.legend)==0) :
          self.legend = a[3]
        if (line.startswith("# selectStr")) :
          self.selectStr = a[3]
        if line.startswith("# avgchan") :
          self.avgchan = a[3]
        if (len(a) > 9) and (line.startswith("C")) :		# new style Lk table, includes all antennas
          ant = int( a[0].strip("C") )
          if ant == self.ant :
            f1 = min(float(a[1]),float(a[2]))
            f2 = max(float(a[1]),float(a[2]))
            DR = float(a[3]) + 1j * float(a[4])
            DL = float(a[5]) + 1j * float(a[6])
            Qpercent = float(a[7])
            Upercent = float(a[8])
            solin.append( [f1,f2,DR,DL,a[9],Qpercent,Upercent] )
        elif (not line.startswith("#")) and (len(a) > 10) :       # old style lk table, one antenna only
            f1 = min(float(a[0]),float(a[1]))
            f2 = max(float(a[0]),float(a[1]))
            DR = float(a[6]) + 1j * float(a[7])
            DL = float(a[8]) + 1j * float(a[9])
            chanfacts = a[12].split(",")
            self.avgchan = chanfacts[3]
            solin.append( [f1,f2,DR,DL,a[12],0.,0.] )
      fin.close()

    # it would be smarter to store values as self.sol = sorted..., but I am temporarily breaking
    #   up everything into f1,f2,DR,DL for compatibility with existing routines
      for s in sorted( solin, key = lambda(x) : x[1] ):	# store solutions in frequency order
        self.f1.append( s[0] )
        self.f2.append( s[1] )    
        self.DR.append( s[2] )
        self.DL.append( s[3] )
        self.lineStr.append( s[4] )
        self.Qpercent.append( s[5] )
        self.Upercent.append( s[6] )
   
  def list(self) :
    return [self.ant, self.file, self.legend]

  def fminmax(self) :
    fmin = self.f1[0]
    fmax = self.f2[-1]
    return [fmin,fmax] 
    
  # --- usually use DRcolor=None, DLcolor=None ---
  def plotComplex( self, p, fstart, fstop, DRcolor='red', DLcolor='blue' ) :
    """plot DR and DL on the complex plane"""
    first = True
    f2prev = 0.
    DRmark = "o"
    DLmark = "o"
    if not DRcolor :
      DRcolor = self.color
      DRmark = "."
    if not DLcolor :
      DLcolor = self.color
      DLmark = "D"
    p.plot(0.,0.,marker="+", color="green", linewidth=3, markersize=10 )
    for f1,f2,DR,DL in zip( self.f1, self.f2, self.DR, self.DL ) :
      msize = 5
      #if (f2-f1) > 0.24 :
      #  msize = 8
      favg = 0.5 * (f1 + f2)
      if (favg > fstart) and (favg < fstop) and (numpy.abs(DR) != 0.) :
        if first :     # plot dot and legend
          if not DRcolor:
            p.plot( DR.real, DR.imag, marker=DRmark, color=DRcolor, markersize=msize, \
              label=self.legend )
          else :	# skip legend if DR (and presumably DL) color is specified
            p.plot( DR.real, DR.imag, marker=DRmark, color=DRcolor, markersize=msize)
          p.plot( DL.real, DL.imag, marker=DLmark, color=DLcolor, markersize=msize )
          first = False
        else :         
          #if f1 == f2prev :			# connect with line
          if True :
            p.plot( [xRprev,DR.real], [yRprev,DR.imag], marker=DRmark, color=DRcolor, \
             markersize=msize, linestyle='solid', linewidth=1 )
            p.plot( [xLprev,DL.real], [yLprev,DL.imag], marker=DLmark, color=DLcolor, \
             markersize=msize, linestyle='solid', linewidth=1 )
          else :                    # don't connect with line  
            p.plot( DR.real, DR.imag, marker=DRmark, color=DRcolor, markersize=msize )
            p.plot( DL.real, DL.imag, marker=DLmark, color=DLcolor, markersize=msize )
        xRprev = DR.real
        yRprev = DR.imag
        xLprev = DL.real
        yLprev = DL.imag
        f2prev = f2
    if DRcolor :
      p.text( .14, .11, "C%d" % self.ant,fontsize=10, horizontalalignment="right" )
    #p.plot(0.,0.,marker="+", color="green", linewidth=2, markersize=10 )

    p.legend( loc=0, prop={'size':6}, numpoints=1 )
         
  def panel(self, p, type, lk, fstart, fstop, ShowLegends=True ) :
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
          if type == "phs" : #or type == "amp" :
            p.plot( f, y, marker='o', color=self.color, markersize=3 )    # make dot
        # draw lines for amp
          else :
            if ((f1 > fstart) and (f1 < fstop)) or ((f2 > fstart) and (f2 < fstop)) :
              if f1 == f2prev :
                if first :
                  print self.legend
                  p.plot( [fprev,f], [yprev,y] , color=self.color, \
                    linestyle='solid', linewidth=1, label=self.legend )
                  first = False
                else :
                  p.plot( [fprev,f], [yprev,y] , color=self.color, \
                    linestyle='solid', linewidth=1 )
            fprev = f
            f2prev = f2
            yprev = y
    if ShowLegends :
      p.legend( loc=0, prop={'size':6} )
    

# ----------------------------------------------------------------------------------------------------- #
# the LkSet object - a collection of multiple leakage objects, for various antennas, frequencies, sources

class LkSet:
  """LkSet object contains one or more Leak objects"""
 
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
      LkSet.loadAll(self, masterListFile=masterList )

  # --- add one Leak to LkSet object ---
  def addLeak(self, file, antenna, legend, color, marker ) :
    newLeak = Leak( file, antenna, legend, color, marker )
    self.LeakList.append( newLeak )
    self.plotList.append( True )

  def loadAll(self, masterListFile ) :
    """Read in multiple leakage files to LkSet object"""
    color = [ "red", "blue", "green", "chartreuse", "orangered", \
              "aqua", "fuchsia", "gray", "lime", "maroon", "navy", \
              "olive", "orange", "silver", "teal", "black" ]
  # --- used these colors for mxrbias leakage figure --- #
    # color = [ "red", "MediumAquamarine", "LightSeaGreen", "Aqua", "DarkCyan", \
    #           "Teal", "DarkTurquoise", "gray", "lime", "maroon", "navy", \
    #           "olive", "orange", "silver", "teal", "black" ]
    marker = [ "o", "D", "v", "^", "s", "h", "d", \
               "o", "D", "v", "^", "s", "h", "d", \
               "o", "D" ]
    fin = open( masterListFile, "r" )
    ncolor = 0
    for line in fin :
      if len(line) > 0 :                    # skip blank lines
        truncate = line.find("#")           # skip anything after "#"
        a = line[0:truncate].split()
        b = []
        if len(a) > 1 :
          b = a[1].split(",")
          print "skip ants: ", b
        if len(a) > 0 :  
          for nant in range(1,16) :
            if a[0][-1] == "*" :            # old style lk files, one per antenna
              filename = a[0][0:-1] + str(nant)
            else :                          # new style Lk file, all antennas in 1 file
              filename = a[0]
            legend = ""
            if len(a) > 1 : legend = a[1]
            if "%d" % nant in b :
              print "SKIPPING ANT %d" % nant
            else :
              LkSet.addLeak(self, filename, nant, legend, color[ncolor], marker[ncolor] )
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

  def ampAll(self, f1=0., f2=0., type="amp", amax=.25, antList=allAnts, outfile="AmpLeaks.pdf" ) :

    ShowPolfits = True
    ShowLegends = True

    pyplot.ioff()
    if f1 == 0. :
      [f1, f2] = LkSet.xlimits( self )          # default is to find freq limits in the data
    print "frequency limits: %.3f - %.3f GHz" % (f1,f2)
    ymin = 0.
    ymax = amax
    pp = PdfPages( outfile )
    if type == "phs" :
      ymin = -180.
      ymax = 180.
      pp = PdfPages( 'PhsLeaks.pdf' )
    pyplot.clf()

    for ant in antList :
      pL = pyplot.subplot(2, 1, 1)    # DL in upper panel
      pL.axis( [f1, f2, ymin, ymax], size=8 )  # adding size does nothing!!
      pL.grid(True)
      pL.set_ylabel('leakage amplitude', fontsize=10)

      if ShowPolfits : 
        f = []
        y = []
        fin = open("/o/plambeck/pol/LkLib/polfits/C%d.polfit" % ant, "r")
        for line in fin :
          a = line.split()
          if not line.startswith("#") :
            f.append( float(a[0]) )
            y.append( float(a[1]) ) 
          else :
            polname = a[1]
        labeltext = None
        #if ShowLegends : 
        #  labeltext = "%s expected" % polname
        pL.plot( f, y, linestyle="dashed", label=labeltext )

      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DL for antenna %d" % ant
          Leak.panel(pL, type, "DL", f1, f2, ShowLegends=ShowLegends )
      pyplot.title("C%d DL" % ant)

      pR = pyplot.subplot(2, 1, 2)    # DR in lower panel
      pR.axis( [f1, f2, ymin, ymax] )
      pR.grid(True)
      pR.set_ylabel('leakage amplitude', fontsize=10)
      pR.set_xlabel('frequency (GHz)', fontsize=10)
      if ShowPolfits :
        labeltext = None
        #if ShowLegends : 
        #  labeltext = "%s expected" % polname
        pR.plot( f, y, linestyle="dashed", label=labeltext )
        pR.legend( loc=0, prop={'size':6} )
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DR for antenna %d" % ant
          Leak.panel(pR, type, "DR", f1, f2, ShowLegends=False )
      pyplot.title("C%d DR" % ant)
      pyplot.savefig( pp, format='pdf', bbox_inches='tight'  )
      pyplot.clf()
      #pyplot.savefig( pp, format='pdf', bbox_inches='tight'  )
    pp.close()

# for complex plot, f1 and f2 could in principle be used to select range of points
#    but this is not implemented yet

  def complexAll(self, f1=0., f2=0., amax=.16, nrows=1, ncols=1, antList=allAnts ) :
    pyplot.ioff()
    pp = PdfPages( 'ComplexLeaks.pdf' )
    scale = 10./math.sqrt(ncols*nrows)
    if f1 == 0. :
      [f1, f2] = LkSet.xlimits( self )          # default is to find freq limits in the data
    print "frequency limits: %.3f - %.3f GHz" % (f1,f2)
    ymin = -1.*amax
    ymax = amax
    npanel = 0
    for ant in antList :
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
      #pyplot.title("C%d DR (circles, solid) and DL (diamonds, dashed)" % ant, fontdict={'fontsize': scale})
      if (npanel == nrows*ncols) or (ant == antList[-1] ) :
        pyplot.savefig( pp, format='pdf' )
    pp.close()

# first file in LeakList is Lktemplate
# for each freq interval in Lktemplate that is within range f1-f2, find all other entries in LeakList
#   which match within 10% of the freq range; append these to list; form avg of list; plot scatter
#   relative to the avg, 15 plots per page; compute rms

  def complexCmp(self, fmin=0., fmax=300., amax=.05, nrows=4, ncols=4, delta=True, antList=allAnts ) :
    Lktemplate = self.LeakList[0]
    DRmark = "."
    DLmark = "x"
    pyplot.ioff()
    pp = PdfPages( 'ComplexCmp.pdf' )
    fout = open( "ComplexCmp.dat", "w" )
    fout2 = open( "ComplexCmpDelta.dat", "w" )
    fout3 = open( "LkAvg", "w" )
    for f1ref,f2ref,lineStr in zip ( Lktemplate.f1, Lktemplate.f2, Lktemplate.lineStr ) :
      if (f1ref > fmin) and (f2ref < fmax) :
        finterval = f2ref-f1ref
        fout.write("\n ******************************************************************* \n")
        fout.write("(%.3f,%.3f)\n" % (f1ref,f2ref))
        print f1ref,f2ref,finterval
        npanel = 0

      # Process one antenna at a time 
        for ant in antList :
          print ""
          fout.write("\n")
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

        # redefine color map so each LeakFile has unique color
        # ... this is a FLAWED procedure if there will be only one legend for all antennas,
        # ... since colors could be different for different antennas
          cmap = pyplot.get_cmap("brg")
          numcolors = len(collist)
          for n in range(0,numcolors) :
            collist[n] = cmap(n/float(numcolors))

        # summarize raw data (without subtracting the mean) in file
          DRmean = numpy.mean(DRlist)
          DLmean = numpy.mean(DLlist)
          rmsDR = numpy.std( DRlist ) 
          rmsDL = numpy.std( DLlist ) 

          for DR,DL,legend in zip( DRlist, DLlist, leglist ) :
            print "%3d  (%+5.3f%+5.3fj) %5.3f   (%+5.3f%+5.3fj) %5.3f   %s" % \
             ( ant, DR.real, DR.imag, abs(DR), DL.real, DL.imag, abs(DL), legend )
            fout.write("%3d  (%+5.3f%+5.3fj) %5.3f   (%+5.3f%+5.3fj) %5.3f   %s\n" % \
             ( ant, DR.real, DR.imag, abs(DR), DL.real, DL.imag, abs(DL), legend ) )
            
            fout2.write(" C%d  %.4f  %.4f   %.4f    deltaDR\n" % 
             (ant, numpy.real(DR-DRmean), numpy.imag(DR-DRmean), abs(DR-DRmean)))
            fout2.write(" C%d  %.4f  %.4f   %.4f    deltaDL\n" % 
             (ant, numpy.real(DL-DLmean), numpy.imag(DL-DLmean), abs(DL-DLmean)))

          print "%3d  (%+5.3f%+5.3fj) %5.3f   (%+5.3f%+5.3fj) %5.3f   %s" % \
           ( ant, DRmean.real, DRmean.imag, rmsDR, DLmean.real, DLmean.imag, rmsDL, "AVG" )
          fout.write( "%3d  (%+5.3f%+5.3fj) %5.3f   (%+5.3f%+5.3fj) %5.3f   %s\n" % \
           ( ant, DRmean.real, DRmean.imag, rmsDR, DLmean.real, DLmean.imag, rmsDL, "AVG" ))
          percentQ = 0.
          percentU = 0.
          RLgainRatio = 0.
          fout3.write("C%02d %8.3f %8.3f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f   %s   %5.3f\n" % \
           ( ant, f1ref, f2ref, DRmean.real, DRmean.imag, DLmean.real, \
             DLmean.imag, percentQ, percentU, lineStr, RLgainRatio) )


        # if delta=True, subtract average before plotting
          if (delta) :
            DRlist = DRlist - DRmean
            DLlist = DLlist - DLmean
          DRref = numpy.mean(DRlist)  # circles will be centered on DRref 
          DLref = numpy.mean(DLlist)
        
        # plot this panel
          npanel= npanel + 1
          p = pyplot.subplot(nrows, ncols, npanel, aspect='equal')    # DL,DR in one panel
          p.tick_params( axis='both', which='major', labelsize=4 )
          #p.set_ticks( [-.1, 0., .1] )
          p.axis( [-1.*amax, amax, -1.*amax, amax] )
          p.plot( [0.,0.], [-amax,amax], "k--", linewidth=0.2 )
          p.plot( [-amax,amax], [0.,0.], "k--", linewidth=0.2 )
          #p.grid(True)

        # Draw circles showing scatter; center on mean position, or on 0,0 if delta==True
          circ1 = Circle( (DRref.real,DRref.imag), rmsDR, linestyle="dotted", color="r", fill=False )
          #p.add_patch(circ1)
          circ2 = Circle( (DLref.real,DLref.imag), rmsDL, linestyle="dotted", color="b", fill=False )
          #p.add_patch(circ2)

          #p.text(.07,.07, "DR = %5.3f%+5.3fj" % (DRmean.real, DRmean.imag), horizontalalignment='left',transform=p.transAxes, size=6, color="r" )
          #p.text(.07,.15, "DL = %5.3f%+5.3fj" % (DLmean.real, DLmean.imag), horizontalalignment='left',transform=p.transAxes, size=6, color="b" )
          #p.text(.07,.17, "DL %5.3f" % rmsDL, horizontalalignment='left',transform=p.transAxes, size=7, color="b" )
          p.text(.93,.85,"C%d" % ant, horizontalalignment='right',transform=p.transAxes, size=7)

          for DR,legend,color in zip( DRlist,leglist,collist ) :
            p.plot( DR.real, DR.imag, marker=DRmark, color=color, markeredgecolor=color, markersize=4, \
              label=legend )
          for DL,color in zip( DLlist,collist) :
            p.plot( DL.real, DL.imag, marker=DLmark, color=color, markeredgecolor=color, markersize=4 )
          if ant == 13 :
            p.legend(bbox_to_anchor=(5.5,0.9), prop={'size':3})
            p.text(5.52,0.92,"(%.3f,%.3f)" % (f1ref,f2ref), horizontalalignment='right',transform=p.transAxes, size=7)
        pyplot.savefig( pp, format='pdf' )
        pyplot.clf()
    pp.close()
    fout.close()
    fout2.close()
    fout3.close()
         

# average together leakages in the LkSet object, write out avg to new LkFile; optionally, add offset
# remember that each leak inside LeakList covers just one antenna

  def newLk( self, newLkFile, DRoffset=0.+0j, DLoffset=0.+0j ) :
    fout = open( newLkFile, "w" )
    fout.write("# input data : %s\n" % "AVG"  )
    percentU = percentQ = 0.
    for lineStr,f1,f2 in zip (self.LeakList[0].lineStr, self.LeakList[0].f1, self.LeakList[0].f2) :
      print lineStr
      fout.write("#\n")
      for ant in range(1,16) :
        DRlist = []
        DLlist = []
        for Lk in self.LeakList :
          if Lk.ant == ant :
            for lineStr1,DR1,DL1 in zip( Lk.lineStr, Lk.DR, Lk.DL ) :
              if (lineStr1 == lineStr) and (abs(DR1) > 0.) and (abs(DL1) > 0.) :
                DRlist.append( DR1 )
                DLlist.append( DL1 )
                print "... ant %d - appending data from %s" % ( ant, Lk.legend )
        DRmean = numpy.mean(DRlist) + DRoffset
        DLmean = numpy.mean(DLlist) + DLoffset
        print ant, DRlist, DRmean, DLlist, DLmean
        fout.write("C%02d %8.3f %8.3f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f    %s\n" % \
            ( ant, f1, f2, DRmean.real, DRmean.imag, DLmean.real, \
            DLmean.imag, percentQ, percentU, lineStr) )
    fout.close()

# average together leakages for each antenna in the LkSet object; 
# offset real and imag parts of the leakage with a Gaussian random error
#   (error is the same for real and imaginary, for all ants)
# apply Miriad constraint?; write new Lk file 

  def GaussLk( self, newLkFile, sigma ) :
    fout = open( newLkFile, "w" )
    fout.write("# input data : %s\n" % "+ Gaussian error"  )
    percentU = percentQ = 0.
    for lineStr,f1,f2 in zip (self.LeakList[0].lineStr, self.LeakList[0].f1, self.LeakList[0].f2) :
      print lineStr
      fout.write("#\n")
      for ant in range(1,16) :
        DRlist = []
        DLlist = []
        for Lk in self.LeakList :
          if Lk.ant == ant :
            for lineStr1,DR1,DL1 in zip( Lk.lineStr, Lk.DR, Lk.DL ) :
              if (lineStr1 == lineStr) and (abs(DR1) > 0.) and (abs(DL1) > 0.) :
                DRlist.append( DR1 )
                DLlist.append( DL1 )
                print "... ant %d - appending data from %s" % ( ant, Lk.legend )
        if len(DRlist) > 0 :
          DRmean = numpy.mean(DRlist)
          DLmean = numpy.mean(DLlist)
          DRnew = random.gauss(numpy.real(DRmean), sigma) + random.gauss(numpy.imag(DRmean), sigma) * 1j
          DLnew = random.gauss(numpy.real(DLmean), sigma) + random.gauss(numpy.imag(DLmean), sigma) * 1j
        else :
          DRnew = 0. + 0j
          DLnew = 0. + 0j
        print ant, DRnew, DLnew
        fout.write("C%02d %8.3f %8.3f %8.3f %6.3f %8.3f %6.3f %8.3f %6.3f    %s\n" % \
            ( ant, f1, f2, DRnew.real, DRnew.imag, DLnew.real, \
            DLnew.imag, percentQ, percentU, lineStr) )
    fout.close()
     
# Figure out how much leakage changes the XYphase calibration
# assume that L and R are in phase
# then compute measured angles of L and R taking into account the leakge
# then compute difference in these phases - it turns out to be very small
# note: this is specialized for DSB data where there is only one DR,DL for each ant

  def XYphsOffset( self ) :
    for Lk in self.LeakList :
      VRmeas = numpy.angle( 1. + Lk.DR[0], deg=True) 
      VLmeas = numpy.angle( 1. + Lk.DL[0], deg=True)
      #print Lk.ant, VRmeas, VLmeas 
      print "ant %d  %.2f  %.2f  XYphsOffset = %.2f" % (Lk.ant, VLmeas, VRmeas, VLmeas-VRmeas)

# ==== special to make JAI figure ==== #
  def ampJAI2(self, f1=0., f2=0., type="amp", amax=.25, antList=allAnts, outfile="AmpLeaks.pdf" ) :
    ShowPolFits = True
    if f1 == 0. :
      [f1, f2] = LkSet.xlimits( self )          # default is to find freq limits in the data
    fig = pyplot.figure( figsize=(10,10), facecolor='white' )
    xstart = .05
    ystart = .05 
    delx = .14
    dely = 0.10
    ymin = 0.0001
    ymax = amax - .0001

    for ant in antList :

      f = []
      y = []
      fin = open("polfits/C%d.polfit" % ant, "r")
      for line in fin :
        a = line.split()
        if not line.startswith("#") :
          f.append( float(a[0]) )
          y.append( float(a[1]) ) 
        else :
          polname = a[1]

      pR = fig.add_axes( [xstart, ystart, delx, dely], autoscale_on=True )
      pR.tick_params( labelsize=6 )
      pR.axis( [f1, f2, ymin, ymax], size=8 )  # adding size does nothing!!
      pR.grid(True)
      pR.text( 257, .2, "C%d DR" % ant, horizontalalignment='right', fontsize=10 )
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DR for antenna %d" % ant
          Leak.panel(pR, type, "DR", f1, f2, ShowLegends=False )

      pL = fig.add_axes( [xstart, ystart+dely, delx, dely], autoscale_on=True )
      pL.tick_params( labelsize=6 )
      pL.axis( [f1, f2, ymin, ymax], size=8 )  # adding size does nothing!!
      pL.grid(True)
      pL.set_xticklabels( [] )
      pL.plot( f, y, linestyle="dashed" ) # , label="%s expected" % polname )
      pL.text( 257, .2, "C%d DL" % ant, horizontalalignment='right', fontsize=10 )
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DL for antenna %d" % ant
          Leak.panel(pL, type, "DL", f1, f2, ShowLegends=False )
      pR.plot( f, y, linestyle="dashed" ) # , label="%s expected" % polname )

      xstart = xstart + .17 
      ystart = ystart + 0.
      if (xstart+delx) > 1. :
        xstart = 0.05
        ystart = ystart + .25

    pyplot.show() 

# ==== special to make JAI figure, plotting R and L on same panel ====  #

  def ampJAI(self, f1=205., f2=259., type="amp", amax=.2, antList=allAnts, outfile="AmpLeaks.pdf" ) :
    ShowPolFits = True
    if f1 == 0. :
      [f1, f2] = LkSet.xlimits( self )          # default is to find freq limits in the data
    fig = pyplot.figure( figsize=(10,10), facecolor='white' )
    xstart = .1
    ystart = .05 + 4.* 0.18
    delx = .25
    dely = 0.15
    ymin = 0.0
    ymax = amax 

    for ant in antList :

      f = []
      y = []
      fin = open("polfits/C%d.polfit" % ant, "r")
      for line in fin :
        a = line.split()
        if not line.startswith("#") :
          f.append( float(a[0]) )
          y.append( float(a[1]) ) 
        else :
          polname = a[1]

      pR = fig.add_axes( [xstart, ystart, delx, dely], autoscale_on=True )
      pR.tick_params( labelsize=10, width=0.5 )
      pR.axis( [f1, f2, ymin, ymax] )  # adding size does nothing!!
      pR.grid(True)
      pR.text( 254.5, .175, "C%d" % ant, horizontalalignment='center', verticalalignment='center', fontsize=12 )
      for Leak in self.LeakList :
        if Leak.ant == ant :
          print "plotting DR for antenna %d" % ant
          Leak.color = 'red'
          Leak.panel(pR, type, "DR", f1, f2, ShowLegends=False )
          print "plotting DL for antenna %d" % ant
          Leak.color = 'blue'
          Leak.panel(pR, type, "DL", f1, f2, ShowLegends=False )
      pR.plot( f, y, linestyle="dashed", color='black' ) # , label="%s expected" % polname )

    # this is a lame way of adding labels, but I couldn't manage to turn off the dumb frame when
    #   I generated a separate axes to cover the entire figure and labeled it
      if ant == 14 :
        pR.set_xlabel( "frequency (GHz)", fontsize=10)
      if ant == 7 :
        pR.set_ylabel( "leakage amplitude", fontsize=10)


      xstart = xstart + .3 
      ystart = ystart + 0.
      if (xstart+delx) > 1. :
        xstart = 0.1
        ystart = ystart - .18

    pyplot.show() 

# ----------------------------------------------------------------------------------------------------- #
# wrapper routines - all deal with a list of leakages specifed as LkList

def plotAmps( LkList, f1=0., f2=0., amax=0.25, antList=allAnts, outfile="AmpLeaks.pdf" ) :
  p = LkSet( LkList )
  p.ampAll( f1=f1, f2=f2, amax=amax, antList=antList, outfile=outfile ) 

def plotPhases( LkList, f1=0., f2=0., antList=allAnts ) :
  p = LkSet( LkList )
  p.ampAll( f1=f1, f2=f2, type='phs', antList=antList ) 

def plotComplex( LkList, f1=0., f2=0., antList=allAnts ) :
  p = LkSet( LkList )
  p.complexAll( f1=f1, f2=f2, amax=.16, nrows=1, ncols=1, antList=antList ) 

def cmpComplex( LkList, f1=0., f2=300. ) :
  p = LkSet( LkList )
  p.complexCmp( fmin=f1, fmax=f2, amax=.05, nrows=4, ncols=4) 

def makeSS( LkFile, SSfile ) :
  oneLeak = Leak( LkFile, 1, "", "black", "o")    # dummy color and markers are required, sadly
  ss = RM.SS()	                                  # empty Stokes Spectrum
  ss.LkFile = LkFile
  ss.selectStr = oneLeak.selectStr 
  ss.visFile = "from LkFile"
  ss.UT = 0.
  ss.parang = 0.
  ss.HA = 0.
  nfreqs = len(oneLeak.f1)
  ss.strList = oneLeak.lineStr
  ss.f1 = numpy.array( oneLeak.f1 )
  ss.f2 = numpy.array( oneLeak.f2 )
  ss.I = 100.*numpy.ones( nfreqs )
  ss.rmsI = numpy.ones( nfreqs )
  ss.Q = numpy.array( oneLeak.Qpercent )
  ss.rmsQ = 0.1 * numpy.ones( nfreqs )
  ss.U = numpy.array( oneLeak.Upercent )
  ss.rmsU = 0.1 * numpy.ones( nfreqs )
  ss.V = numpy.zeros( nfreqs )
  ss.rmsV = numpy.zeros( nfreqs )
  ss.fitPARM( 225. )
  ss.plot()
#  ss.dump( None )
  fout = open( SSfile, "ab" )
  pickle.dump( ss, fout )
  fout.close()

# opens file with single column of numbers, computes std of it
def deltaRMS( deltaFile, column=0 ) :
  delta = []
  fin = open( deltaFile, "r" )
  for line in fin :
    a = line.split()
    delta.append( float(a[column] ) )
  fin.close()
  d = numpy.array(delta)
  print "average: ", numpy.average(delta)
  print "std dev: ",  numpy.std(d)
