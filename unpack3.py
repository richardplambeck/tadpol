# unpack3 processes OpticsBench files

import math
import numpy

# unpack converts raw data file in lock-in format to more readable data file
# comment lines begin with '#' and are copied directly to output file
# fields within each line are:
#   0 - time in seconds since ...
#   1 - freq in GHz
#   2 - amplitude in V
#   3 - phase in degrees
#   4 - position in volts (0-4.516 = 0-32768 counts in carma system)
#   5 - LOpower in mV (from detector on directional coupler before mixer)
# this is an intermediate step in creating transmission files
#
def unpack( infile, outfile ) :
    fin = open( infile,"r" )
    fout = open( outfile,"w" )
    First = True
    for line in fin :
      if line.startswith("#") :      # copy comment line 
        fout.write("%s" % line)
      else :
        a = line.split(",")
        freq = float(a[1])
        amp = float(a[2]) 
        phs = float(a[3]) 
        pos = float(a[4]) 
        LOpwr = float(a[5]) 
        time = float(a[0])
        if First :
          t0 = time
          First = False
        fout.write("%8.2f  %8.6f  %8.3f  %8.5f  pos %5.3f  %8.1f \n" % (freq, amp, phs, LOpwr, pos, time-t0) )
    fin.close()
    fout.write("# end\n")    # needed by trans, below
    fout.close()

# interpret data file created by unpack
# refpos = position of reference, usually hole in puck transporter (but could be alumina puck?)
# Ppos = position of puck 
# currently, blocks of data at given frequency must be separated by comment lines
# set minline=1 for 4pucksH,K,L,M.raw
#
def trans( infile, outfile, refpos, Ppos, tol=0.2, minline=0 ) :
    fin = open( infile,"r")
    print "\n%s" % outfile
    fout = open( outfile,"w")
    fout2 = open( "repeat.dat", "w" )
    fout2.write("# %s\n" % infile )
    fout3 = open( "refamps", "w" )
    nline = 0
    nref = 0
    nP = 0
    npwr = 0
    ampref = 0.
    phsref = 0.
    ampP = 0.
    phsP = 0.
    LOpwr = 0.
    Bad = False
    f = []
    t = []
    p = []
    lo = []

    frefsave = []
    amprefsave = []
    phsrefsave = []

    for line in fin :

    # copy all comment lines to output file
      if line.startswith("#") and not line.startswith("# freq:") :   
        fout.write("%s" % line )
        fout3.write("%s" % line )

    # new frequency; write out previously accumulated daa
      if line.startswith("#") :          
        if (nref == 0) or (nP == 0) :    # first line in file, or perhaps some kind of screwup
          continue 
        elif Bad :
          print "bad data - skipping frequency %.2f" % freq
        else :
          frefsave.append( freq )
          trans = (ampP/nP)/(ampref/nref)
          amprefsave.append( ampref/nref ) 
          try :
            transTrue = pow(10.,1.02*math.log10(trans))    # correct for nonlinearity, as measured with Saturation.txt
          except :
            print "MATH ERROR: trans=%.3f" % trans
            transTrue = 0.
          delphs = phsref/nref - phsP/nP
          phsrefsave.append( phsref/nref )
          while (delphs > 180.) :
            delphs = delphs - 360.
          while (delphs < -180.) :
            delphs = delphs + 360.
          f.append( freq )
          t.append( transTrue )
          p.append( delphs ) 
          lo.append( LOpwr/npwr )
          #fout.write("%8.2f  %8.6f  %8.3f  %8.6f\n" % (freq,trans,delphs,LOpwr/npwr))
          
      # rezero all the buffers
        nref = 0
        nP = 0
        npwr = 0
        ampref = 0.
        phsref = 0.
        ampP = 0.
        phsP = 0.
        LOpwr = 0.
        nline = 0  # nline is here just because of the bug that affected 4pucksK-4pucksM
        Bad = False

    # process data line
      else :
        nline = nline+1
        a = line.split()
        freq = float(a[0])

      # if even one entry for this freq has amp < .0001, discard the entire freq
        if float(a[1]) < .0001 :      
          Bad = True
        pos = float(a[5])

      # create avg LOpwr at this freq, for ref or puck lines
        LOpwr = LOpwr + float(a[3]) 
        npwr = npwr + 1

      # process reference line
        if abs(pos-refpos) < tol and (nline > minline) :
          nref += 1

        # deal with case where ref phase jitters between +/-180
          phstmp = float(a[2])
          if nref == 2 and phstmp > phsref + 180. :
            phstmp = phstmp - 360.
          if nref == 2 and phstmp < phsref - 180. :
            phstmp = phstmp + 360.

        # dump difference between ref values to file repeat.dat, in case it is needed
          if nref == 2 and ampref > .0001 :
            ratio = float(a[1])/ampref
            deltaphs = phstmp - phsref
            fout2.write("%8.2f  %8.6f  %8.3f\n" % (freq,ratio,deltaphs))
            if abs(deltaphs) > 5. :
              Bad =True

          ampref = ampref + float(a[1])
          phsref = phsref + phstmp  

      # process puck line (normally, there is only 1 measurement so we don't check for 180 degree phase change))
        elif abs(pos-Ppos) < tol :  
          ampP = ampP + float(a[1])
          phsP = phsP + float(a[2])
          nP = nP + 1

    fin.close()

  # reorder the data by frequency before writing it out
    iarr = numpy.argsort( f )
    for i in range(0,len(iarr) ) :
        j = iarr[i]
        fout.write("%8.2f  %8.6f  %8.3f  1.0   %8.6f  %8.6f\n" % (f[j],t[j],p[j],lo[j],t[j]*t[j]))
             # the 1.0 is in the uncertainty column
    fout.close()
    fout2.close()

  # now dump out amplitudes of reference, to search for standing waves
  # add angIdeg line and 4th column with 0 so I can read with readTransFile
    iarr = numpy.argsort(frefsave)
    for i in range(0,len(iarr)) :
        j = iarr[i]
        fout3.write("%8.2f  %8.6f  %8.3f  0.0\n" % (frefsave[j],amprefsave[j],phsrefsave[j]) )
    fout3.close()


# average data into chunks, compute rms of each chunk
def binit( infile, frqinterval=0.199, minpts=3 ) :
    fin = open(infile+".dat","r")
    fout = open(infile+".avg.dat","w")
    First = True       # will pick out first data line; write header above this
    fout.write( "# binned data from %s\n" % (infile+".dat") )
    f0 = 0.
    npts = 0
    frq = []
    vec = []
    for line in fin:
      if line.startswith("#") :
        fout.write( "%s" % line )
      else :

        if First :
          fout.write( "# ------------------------------------------------------------------ \n")
          fout.write("#  frq     amp_avg   phs_avg    vec_std     pwr_avg   pwr_std   navg\n")
          First = False

        a = line.split()
        if (float(a[0]) - f0) >= frqinterval :
          if len(vec) >= minpts :
            vecavg = numpy.average( vec )
            if len(vec) > 1 :
              vecstd = numpy.std(vec,ddof=1)/math.sqrt(float(npts))
              pwrstd = numpy.std(pwr,ddof=1)/math.sqrt(float(npts))
            else :
              vecstd = 1.
              pwrstd = 1.
            for i in range(0,len(vec)) :
              if abs(vec[i] - vecavg) > 2.*math.sqrt(npts)*vecstd :
                print "noisy point at %.2f" % frq[i], abs(vec[i] - vecavg), 3.*math.sqrt(npts)*vecstd
                #vec.remove(i)
              #npts = len(vec)
            frqavg = numpy.average( frq )
            fout.write("%8.2f  %8.6f  %8.3f  %10.6f    %8.6f  %8.6f   %2d\n" % (frqavg, numpy.absolute(vecavg), \
              numpy.angle(vecavg,deg=True), vecstd, numpy.average(pwr), pwrstd, npts)) 
          else :
            fout.write("# missing bin; npts = %d, too few to compute uncertainty\n" % npts)
          f0 = float(a[0])
          npts = 0
          frq = [] 
          vec = []
          pwr = []
        frq.append( float(a[0]) )
        vec.append( float(a[1]) * numpy.exp( 1j*math.radians(float(a[2])) ) )   # convert to real+j*imag
        pwr.append( pow( float(a[1]), 2) )
        npts = npts + 1

  # remember to print out the last point
    if npts >= minpts :
      vecavg = numpy.average( vec )
      frqavg = numpy.average( frq )
      if npts > 1 :
        vecstd = numpy.std(vec,ddof=1)/math.sqrt(float(npts))
        pwrstd = numpy.std(pwr,ddof=1)/math.sqrt(float(npts))
      else :
        vecstd = 1.
        pwrstd = 1.
      fout.write("%8.2f  %8.6f  %8.3f  %10.6f    %8.6f  %8.6f   %2d\n" % (frqavg, numpy.absolute(vecavg), \
        numpy.angle(vecavg,deg=True), vecstd, numpy.average(pwr), pwrstd, npts)) 
    else :
      fout.write("# missing bin; npts = %d, too few to compute uncertainty\n" % npts)
    fin.close()
    fout.close()
