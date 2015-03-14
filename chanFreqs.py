# calculate and plot VLBI frequency setup
# begin with synthesizer freq and harmonic multiples, compute exact LO1 freq
# then figure out exact LO2 freq for each CARMA band, from fsky at MIDPOINT of CARMA band
#   (just as it's done with configastroband)
# once LO1, LO2, nsb are fixed, then specifying any base frequency gives sky freq
# thus, I can compute sky freq at midpoint of each dbe channel
# then compute precise sky freq at midpoint of each DBE channel
# begin with LO1 freq, midpoint of CARMA band, and sideband of LO1 - compute LO2, fblock
#
import math
import numpy
from matplotlib import pyplot
from matplotlib.patches import Rectangle

# --- specify LO1 freq in MHz --- #  
# fLO1MHz = 222100.00000004997
   # this is the value closest to 222.100 GHz; synth = 1175.449735450

#fLO1MHz = 222099.999915   
   # this is the value that gives closest integer Hz; synth = 1175.449735000

fLO1MHz = 221099.9999700000
   # this is closest integer Hz freq to 221.1 GHz (for 5-7 filter);
   # synth = 1170.1587300000

# --- specify center freq of lowest CARMA band --- #
#win1center = 4245
   # note that 4242.5, 4241.25 etc will work also

win1center=5320.

# --- specify CARMA band centers in MHz --- #
#wincenters = [ win1center, win1center+480., win1center+480.+590., win1center+480.+590.+480. ] 
  # case A: best fits DBE chans into the available bandwidth
  # however, chans above and below 5 GHz are on different 32 MHz grids
#wincenters = [ win1center, win1center+480., win1center+480.+640., win1center+480.+640.+480. ] 
  # case B: keeps all DBE chans on one 32 MHz ladder
#wincenters = [ win1center, win1center+480., win1center+2*480.+100., win1center+3*480.+100. ] 
  # case C: this is the correct answer for 4-6 GHz IF ! 
wincenters = [ win1center, win1center+480., win1center+2*480., win1center+3*480. ] 
  # case D: use for 5-7 GHz

def computeLO1( fsynthMHz=1175.449735450, nharm=9, mharm=7 ) :
#def computeLO1( fsynthMHz=1175.449735000, nharm=9, mharm=7 ) :
  xfreq = nharm * fsynthMHz - 10.
  gunnfreq = mharm * xfreq + 50.
  LO1 = 3.*gunnfreq
  print "LO1 calculation:"
  print "... synth = %.10f MHz" % fsynthMHz
  print "... xband = %.10f MHz" % xfreq
  print "... gunn = %.10f MHz" % gunnfreq
  print "... LO1 = %.10f MHz" % LO1
  return LO1

# LO2 must be a multiple of (10 MHz)/2^{18} 
def exactLO2( targetLO2MHz ) :
  multiplier = round(targetLO2MHz*1.e6/(1.e7/pow(2,18)))
  LO2actual = (multiplier * 1.e7/pow(2,18))/1.e6
  print "... LO2 target = %.9f, actual = %.11f MHz" % (targetLO2MHz, LO2actual )
  return LO2actual

def doit() :
  findLO2( fLO1MHz+4230. )

# given MIDBAND sky freq (as used at CARMA)  and LO1 frequencies, return fLO2, etc. 
def findLO2( fskyMHzMidband, fLO1=fLO1MHz ) :
  print " "
  print "... fsky = %.9f MHz" % fskyMHzMidband
  print "... fLO1 = %.9f MHz" % fLO1
  if fskyMHzMidband > fLO1 : nsb = 1
  else : nsb = -1
  print "... nsb = %d" % nsb
  fIF = nsb * (fskyMHzMidband - fLO1)
  print "... fIF = %.3f" % fIF
  if fIF < 5000. : 
    fIF2 = fIF
    fblock = 0.
  else : 
    fIF2 = 10000. - fIF
    fblock = 10000.
  print "... fIF2 = %.3f" % fIF2
  print "... fblock = %.1f MHz" % fblock
# NOTE: this should be 3500, but I want to use override USB option for IF=3250.
  #if fIF2 > 3500. :                               
  if fIF2 > 3000. :                               
    fLO2 = fIF2 - 750.
    nsb2 = 1   # postiive nsb2 means fIF2 = fLO2 + base
  else : 
    fLO2 = fIF2 + 750.
    nsb2 = -1  # negative nsb2 means fIF2 = fLO2 - base
  print "... fLO2 = %.3f" % fLO2
  LO2actual = exactLO2(fLO2)
    # LO2 must be a multiple of (10 MHz)/2^{18} 
  print "... fLO2actual = %.3f MHz" % LO2actual
  print "... nsb2 = %d" % nsb2
  return nsb, fblock, fLO2, nsb2

# compute sky freq from baseband freq, for specified fLO2, etc.
def fsky( fbase, fLO2, nsb2, fblock, fLO1, nsb ) :
  if nsb2 > 0 : fIF2 = fLO2 + fbase
  else : fIF2 = fLO2 - fbase
  if fblock > 1. : fIF = fblock - fIF2
  else : fIF = fIF2
  fsky = fLO1 + nsb * fIF
  #print " ... fbase = %.2f, fIF2 = %.4f, fIF = %.4f, fsky = %.4f" % (fbase, fIF2, fIF, fsky)
  return fsky

def fDBE( nDBE ) :
  return (1024. - nDBE * 32.)

# =========== BEGIN PLOTTING ============ #
# make plot of band layout, table of start,stop freqs for all DBE channels
def makefplot( fLO1=fLO1MHz, wincenters=wincenters, outfile="DBEfreqs.txt"  ) :

  #bandcolors = [ 'plum', 'coral', 'cyan', 'chartreuse' ]
  bandcolors = [ 'coral', 'chartreuse', 'cyan', 'plum' ]
  #fout = open( outfile, "w" )
  #fout.write("CARMA DBE frequencies for:\n")
  #fout.write("  LO1 = " + "{:,.3f}".format(fLO1*1.e6) + " Hz\n" )
  #for [i,wincenter] in enumerate( wincenters ) :
  #  fout.write("  band%d = " % (2*i-1) + "{:,.3f}".format(wincenter*1.e6) + " Hz\n")
  #fout.write("\n\n")

  fig = pyplot.figure()
  p = fig.add_axes( [0,0,1,1], frameon=False, axis_bgcolor='white' )
  p.get_xaxis().set_visible(False)
  p.get_yaxis().set_visible(False)
  p.set_xlim( [0.1,0.9] )
  p.set_ylim( [226000.,228500.] )
  arrowhead = 50.  #length of arrowhead in MHz 

  # outer rectangle indicates target range
  chanblock = Rectangle( (0.30, 226100.), .28, 2000., fill=False, linestyle='dashed')   # target range
  p.add_patch(chanblock)

  # plot forbidden zone
  fsky0 = fLO1+5000.
  fsky1 = fLO1+5050.
  chanblock = Rectangle( (0.30, fsky0), .28, (fsky1-fsky0), \
    fill=False, facecolor='red', hatch='/', edgecolor='none', alpha=0.5)   # plot the CARMA 500 MHz wide band
  p.add_patch(chanblock)

  # annotate key sky freqs
  for fann in [ 226100., fLO1+5000., fLO1+5050., 228100.] :
    p.annotate( "%.0f" % fann, xy=(.3,fann), xytext=(.26,fann), \
      horizontalalignment='right', verticalalignment='center', \
      arrowprops=dict(facecolor='black', width=0.0, headwidth=0, shrink=0.1) )
      #arrowprops=dict(facecolor='black', width=0.3, headwidth=5, shrink=0.08) )

  # go through CARMA bands; plot beamformer channels
  offset = 0.
  rwidth = .04
  nband = -1 

  for fcenter,bandcolor in zip( wincenters, bandcolors)  :
    nband = nband + 2    # this is just the CARMA label for the band (1,3,5,7)
    [nsb, fblock, fLO2, nsb2] = findLO2( fLO1+fcenter, fLO1 )
    fsky0 = fsky( 500., fLO2, nsb2, fblock, fLO1, nsb )
    fsky1 = fsky( 1000., fLO2, nsb2, fblock, fLO1, nsb )

    # plot the CARMA 500 MHz band; arrow shows direction of increasing freq in final IF block that is sent to beamformer
    lblock = 0.32 + offset
    chanblock = Rectangle( (lblock, fsky0), rwidth, (fsky1-fsky0),\
       edgecolor='none',  fill=True, facecolor=bandcolor, alpha=0.4)  
    p.add_patch(chanblock)
    chanblock = Rectangle( (lblock, fsky0), rwidth, (fsky1-fsky0),\
       edgecolor='black',  fill=False, facecolor=bandcolor)  
    p.add_patch(chanblock)
    arrowstart = fsky0 + 0.15 * (fsky1-fsky0)
    if fsky1 > fsky0 : length = 300.
    else: length = -300.
    p.arrow( lblock+rwidth/2., arrowstart, 0., length, lw=1, \
       head_width=0.01, head_length=arrowhead )

    # annotate IF freq
    p.annotate( "IF %.0f" % fcenter, xy=(lblock+rwidth,(fsky0+fsky1)/2.), \
      xytext=((lblock+rwidth+.03),(fsky0+fsky1)/2.), \
      horizontalalignment='left', verticalalignment='center', \
      arrowprops=dict(facecolor='black', width=0.3, headwidth=5, shrink=0.15) )

    # label the band
    p.annotate( "band%d" % nband, xy=(.29,(fsky0+fsky1)/2.), horizontalalignment='right', \
      verticalalignment='center' )

    # now plot the corresponding DBE channels
    lblock = .52 + offset
    if fcenter < 5000. :
      nstart = 15
      nstop = 0
      nstep = -1
    else :
      nstart = 1
      nstop = 16
      nstep = 1
    for nDBE in range( nstart,nstop,nstep ) :
      fsky0 = fsky( fDBE(nDBE)-16., fLO2, nsb2, fblock, fLO1, nsb )
      fsky1 = fsky( fDBE(nDBE)+16., fLO2, nsb2, fblock, fLO1, nsb )
      #fout.write( " -------- " + "{:,.3f}".format(fsky0*1.e6) + "\n" )
      #fout.write( " DBE%d  \n " % nDBE )
      if nDBE != 1 :
        fillTF = True
        fillColor = bandcolor 
      else: 
        fillTF = False
      chanblock = Rectangle( (lblock, fsky0), rwidth, (fsky1-fsky0) , fill=fillTF, \
          facecolor=bandcolor, alpha=0.5  )   # plot the beamformer bands
      p.add_patch(chanblock)

      p.annotate( "%d"% nDBE, xy=( lblock+rwidth/2., (fsky0+fsky1)/2. ), \
        horizontalalignment='center', verticalalignment='center', fontsize='x-small' )
      #if nDBE == 1 :
      #  p.annotate( "DBE 1", xy=(lblock+rwidth/2., (fsky0+fsky1)/2.), \
      #    horizontalalignment='center', verticalalignment='center', fontsize='x-small' )
      #if nDBE == 15 :
      #  p.annotate( "DBE 15", xy=(lblock+rwidth/2., (fsky0+fsky1)/2.), \
      #    horizontalalignment='center', verticalalignment='center', fontsize='x-small' )

      # label sky frequencies at extremes of DBE channel array
      fplot = fsky0
      #if fblock == 10000. : fplot = fsky1
      if nDBE == 15 :
        p.annotate( "%.0f" % fplot, xy=(lblock+rwidth, fplot), xytext=(.61,fplot), \
          horizontalalignment='left', verticalalignment='center', \
          arrowprops=dict(facecolor='black', width=0, headwidth=0, shrink=0.1) )
       
      fplot = fsky1
      #if fblock == 10000. : fplot = fsky0
      if nDBE == 1 :
        p.annotate( "%.0f" % fplot, xy=(lblock+rwidth, fplot), xytext=(.61,fplot), \
          horizontalalignment='left', verticalalignment='center', \
          arrowprops=dict(facecolor='black', width=0, headwidth=0, shrink=0.08) )

    offset = offset + 0.00

    p.annotate( "CARMA correlator", xy=(.34,228300), horizontalalignment='center', \
      verticalalignment='center' )
    p.annotate( "DBE chans", xy=(.54,228300.), horizontalalignment='center', \
      verticalalignment='center' )
    #p.annotate( "sky freqs for", xy=(.68,228320.), horizontalalignment='center', \
    #  verticalalignment='center', fontsize=8  )
    #p.annotate( "LO1 = %.0f" % fLO1, xy=(.68,228270.), horizontalalignment='center', \
    #  verticalalignment='center', fontsize=8  )
  #fout.close()
  pyplot.show()     
