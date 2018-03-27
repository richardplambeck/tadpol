# unpack2 processes OpticsBench files

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
        fout.write("%8.2f  %8.6f  %8.3f  %8.5f  pos %7.3f  %8.1f \n" % (freq, amp, phs, LOpwr, pos, time-t0) )
    fin.close()
    fout.write("# end\n")    # needed by trans, below
    fout.close()

# interpret data file created by unpack
# refpos = position of reference, usually hole in puck transporter (but could be alumina puck?)
# Ppos = position of puck 
# currently, blocks of data at given frequency must be separated by comment lines
# fitFile contains fitting information (number of layers, thickness limits, etc) for this puck;
#    it does not change from run to run, will be copied into output file
# set minline=1 for 4pucksH,K,L,M.raw
#
def trans( infile, outfile, refpos, Ppos, fitFile=None, tol=0.2, minline=0 ) :
    fin = open( infile,"r")
    fout = open( outfile,"w")
    fout2 = open( "repeat.dat", "w" )
    fout2.write("# %s\n" % infile )
    if fitFile :
      fin2 = open( fitFile, "r")
      for line in fin2 :
        fout.write( "%s" % line )
      fin2.close()
    nline = 0
    nref = 0
    nP = 0
    npwr = 0
    ampref = 0.
    phsref = 0.
    ampP = 0.
    phsP = 0.
    LOpwr = 0.
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

    # new frequency; write out previously accumulated daa
      if line.startswith("#") :          
        if (nref == 0) or (nP == 0) :    # first line in file, or perhaps some kind of screwup
          continue 
        elif ampref < .001 :  
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
          print freq, nref, nP, ampref, ampP
          f.append( freq )
          t.append( transTrue )
          p.append( delphs ) 
          lo.append( LOpwr/npwr )
          # fout.write("%8.2f  %8.6f  %8.3f  %8.6f\n" % (freq,trans,delphs,LOpwr/npwr))
          

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

    # process data line
      else :
        nline = nline+1
        a = line.split()
        freq = float(a[0])
        pos = float(a[5])
        LOpwr = LOpwr + float(a[3]) 
        npwr = npwr + 1

      # process reference line
        if abs(pos-refpos) < tol and (nline > minline) :
          nref += 1

        # dump difference between ref values to file repeat.dat, in case it is needed
          if nref == 2 and ampref > .001 :
            ratio = float(a[1])/ampref
            deltaphs = float(a[2]) - phsref
            fout2.write("%8.2f  %8.6f  %8.3f\n" % (freq,ratio,deltaphs))

          ampref = ampref + float(a[1])

        # deal with case where ref phase jitters between +/-180
          phstmp = float(a[2])
          if nref > 1 and phstmp > phsref + 180. :
            phstmp = phstmp - 360.
          if nref > 1 and phstmp < phsref - 180. :
            phstmp = phstmp + 360.
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
    fout3 = open( "refamps", "w" )
    iarr = numpy.argsort(frefsave)
    for i in range(0,len(iarr)) :
        j = iarr[i]
        fout3.write("%8.2f  %8.6f  %8.3f\n" % (frefsave[j],amprefsave[j],phsrefsave[j]) )
    fout3.close()


# average data into chunks, compute rms of each chunk
# also: dump out power to another file
def binit( infile, frqinterval=0.199, minpts=3 ) :
    fin = open(infile+".dat","r")
    fout = open(infile+".avg.dat","w")
    fout2 = open(infile+".pwr.dat","w")
    f0 = 0.
    npts = 0
    frq = 0. 
    vec = []
    for line in fin:
      if line.startswith("#") :
        fout.write( "%s" % line )
      else :
        a = line.split()
        if (float(a[0]) - f0) > frqinterval :
          if npts >= minpts :
            vecavg = numpy.average( vec )
            if npts > 1 :
              vecstd = numpy.std(vec,ddof=1)/math.sqrt(float(npts))
              pwrstd = numpy.std(pwr,ddof=1)/math.sqrt(float(npts))
            else :
              vecstd = 1.
              pwrstd = 1.
            fout.write("%8.2f  %8.6f  %8.3f  %8.3f  %2d\n" % (frq/npts, numpy.absolute(vecavg), \
              numpy.angle(vecavg,deg=True), vecstd, npts)) 
            fout2.write("%8.2f  %8.6f  %8.6f\n" % (frq/npts, numpy.average( pwr ), pwrstd ))
          else :
            fout.write("# npts = %d\n" % npts)
            fout2.write("# npts = %d\n" % npts)
          f0 = float(a[0])
          npts = 0
          frq = 0. 
          vec = []
          pwr = []
        frq = frq + float(a[0])
        vec.append( float(a[1]) * numpy.exp( 1j*math.radians(float(a[2])) ) )   # convert to real+j*imag
        pwr.append( pow( float(a[1]), 2) )
        npts = npts + 1

  # remember to print out the last point
    if npts >= minpts :
      vecavg = numpy.average( vec )
      if npts > 1 :
        vecstd = numpy.std(vec,ddof=1)/math.sqrt(float(npts))
        pwrstd = numpy.std(pwr,ddof=1)/math.sqrt(float(npts))
      else :
        vecstd = 1.
        pwrstd = 1.
      fout.write("%8.2f  %8.6f  %8.3f  %8.3f  %2d\n" % (frq/npts, numpy.absolute(vecavg), \
        numpy.angle(vecavg,deg=True), vecstd, npts)) 
      fout2.write("%8.2f  %8.6f  %8.6f\n" % (frq/npts, numpy.average( pwr ), pwrstd ))
    else :
      fout.write("# npts = %d\n" % npts)
      fout2.write("# npts = %d\n" % npts)
    fin.close()
    fout.close()
    fout2.close()

def leftover1() :

    unpack( "5pucksD.raw", "5pucksD.dat" )
    trans( "5pucksD.dat", "1801_D.dat", 0.39, 1.22, fitFile="1801.fitFile" ) 
    trans( "5pucksD.dat", "1802_D.dat", 0.39, 2.03, fitFile="1802.fitFile" ) 
    trans( "5pucksD.dat", "1803_D.dat", 0.39, 2.82, fitFile="1803.fitFile" ) 
    trans( "5pucksD.dat", "1806_D.dat", 0.39, 3.64, fitFile="1806.fitFile" ) 
    trans( "5pucksD.dat", "rex_D.dat", 0.39, 4.45, fitFile="1805.fitFile" ) 
    #unpack( "5pucksE.raw", "5pucksE.dat" )
    trans( "5pucksE.dat", "1801_E.dat", 0.39, 1.22, fitFile="1801.fitFile" ) 
    trans( "5pucksE.dat", "1802_E.dat", 0.39, 2.03, fitFile="1802.fitFile" ) 
    trans( "5pucksE.dat", "1803_E.dat", 0.39, 2.82, fitFile="1803.fitFile" ) 
    trans( "5pucksE.dat", "1806_E.dat", 0.39, 3.64, fitFile="1806.fitFile" ) 
    trans( "5pucksE.dat", "rex_E.dat", 0.39, 4.45, fitFile="1805.fitFile" ) 

def leftover() :
    #unpack( "5pucksB.raw", "5pucksB.dat" )
    trans( "5pucksB.dat", "1806_B.dat", 0.39, 1.22, fitFile="1806.fitFile" ) 
    trans( "5pucksB.dat", "1807_B.dat", 0.39, 2.03, fitFile="1807.fitFile" ) 
    trans( "5pucksB.dat", "1808_B.dat", 0.39, 2.82, fitFile="1808.fitFile" ) 
    trans( "5pucksB.dat", "1809_B.dat", 0.39, 3.64, fitFile="1809.fitFile" ) 
    trans( "5pucksB.dat", "1810_B.dat", 0.39, 4.45, fitFile="1810.fitFile" ) 

    #unpack( "5pucksC.raw", "5pucksC.dat" )
    trans( "5pucksC.dat", "1804_C.dat", 0.39, 1.22, fitFile="1804.fitFile" ) 
    trans( "5pucksC.dat", "1811_C.dat", 0.39, 2.03, fitFile="1811.fitFile" ) 
    trans( "5pucksC.dat", "1812_C.dat", 0.39, 2.82, fitFile="1812.fitFile" ) 
    trans( "5pucksC.dat", "1813_C.dat", 0.39, 3.64, fitFile="1813.fitFile" ) 
    trans( "5pucksC.dat", "1814_C.dat", 0.39, 4.45, fitFile="1814.fitFile" ) 

    binit( "1801_E" )
    binit( "1802_E" )
    binit( "1803_E" )
    binit( "1806_E" )
    binit( "rex_E" )

    unpack( "5pucksF.raw", "5pucksF.dat" )
    trans( "5pucksF.dat", "1807_F.dat", 0.39, 1.22, fitFile="1807.fitFile" ) 
    trans( "5pucksF.dat", "1808_F.dat", 0.39, 2.03, fitFile="1808.fitFile" ) 
    trans( "5pucksF.dat", "1809_F.dat", 0.39, 2.82, fitFile="1809.fitFile" ) 
    trans( "5pucksF.dat", "1810_F.dat", 0.39, 3.64, fitFile="1810.fitFile" ) 
    trans( "5pucksF.dat", "1812_F.dat", 0.39, 4.45, fitFile="1812.fitFile" ) 
    binit( "1807_F" )
    binit( "1809_F" )
    binit( "1809_F" )
    binit( "1810_F" )
    binit( "1812_F" )

    unpack( "2pucksG.raw", "2pucksG.dat" )
    trans( "2pucksG.dat", "1809_G.dat", 0.39, 1.22, fitFile="1809.fitFile" ) 
    trans( "2pucksG.dat", "rex_G.dat", 0.39, 2.03, fitFile="1805.fitFile" ) 

    unpack( "2pucksH.raw", "2pucksH.dat" )
    trans( "2pucksH.dat", "1809_H.dat", 0.39, 1.22, fitFile="1809.fitFile" ) 
    trans( "2pucksH.dat", "rex_H.dat", 0.39, 2.03, fitFile="1805.fitFile" ) 

    unpack( "2pucksJ.raw", "2pucksJ.dat" )
    trans( "2pucksJ.dat", "1809_J.dat", 0.39, 1.22, fitFile="1809.fitFile" ) 
    trans( "2pucksJ.dat", "rex_J.dat", 0.39, 2.03, fitFile="1805.fitFile" ) 
    binit( "1809_J" )
    binit( "rex_J" )
    unpack( "5pucksK.raw", "5pucksK.dat" )
    trans( "5pucksK.dat", "1801_K.dat", 0.39, 1.22, fitFile="1801.fitFile" ) 
    trans( "5pucksK.dat", "1802_K.dat", 0.39, 2.03, fitFile="1802.fitFile" ) 
    trans( "5pucksK.dat", "1803_K.dat", 0.39, 2.82, fitFile="1803.fitFile" ) 
    trans( "5pucksK.dat", "1806_K.dat", 0.39, 3.64, fitFile="1806.fitFile" ) 
    trans( "5pucksK.dat", "1807_K.dat", 0.39, 4.45, fitFile="1807.fitFile" ) 
    binit( "1801_K" )
    binit( "1802_K" )
    binit( "1803_K" )
    binit( "1806_K" )
    binit( "1807_K" )
    unpack( "5pucksL.raw", "5pucksL.dat" )
    trans( "5pucksL.dat", "1815_L.dat", 0.39, 1.22, fitFile="1815.fitFile" ) 
    trans( "5pucksL.dat", "1816_L.dat", 0.39, 2.03, fitFile="1816.fitFile" ) 
    trans( "5pucksL.dat", "1817_L.dat", 0.39, 2.82, fitFile="1817.fitFile" ) 
    trans( "5pucksL.dat", "1808_L.dat", 0.39, 3.64, fitFile="1808.fitFile" ) 
    trans( "5pucksL.dat", "1810_L.dat", 0.39, 4.45, fitFile="1810.fitFile" ) 
    binit( "1815_L" )
    binit( "1816_L" )
    binit( "1817_L" )
    binit( "1808_L" )
    binit( "1810_L" )
    unpack( "5pucksM.raw", "5pucksM.dat" )
    trans( "5pucksM.dat", "1811_M.dat", 0.39, 1.22, fitFile="1811.fitFile" ) 
    trans( "5pucksM.dat", "1812_M.dat", 0.39, 2.03, fitFile="1812.fitFile" ) 
    trans( "5pucksM.dat", "1813_M.dat", 0.39, 2.82, fitFile="1813.fitFile" ) 
    trans( "5pucksM.dat", "1814_M.dat", 0.39, 3.64, fitFile="1814.fitFile" ) 
    trans( "5pucksM.dat", "1804_M.dat", 0.39, 4.45, fitFile="1804.fitFile" ) 
    binit( "1811_M" )
    binit( "1812_M" )
    binit( "1813_M" )
    binit( "1814_M" )
    binit( "1804_M" )
    unpack( "5pucksN.raw", "5pucksN.dat" )
    trans( "5pucksN.dat", "1818_N.dat", 0.39, 1.22, fitFile="1818.fitFile" ) 
    trans( "5pucksN.dat", "1822_N.dat", 0.39, 2.03, fitFile="1822.fitFile" ) 
    trans( "5pucksN.dat", "1823_N.dat", 0.39, 2.82, fitFile="1823.fitFile" ) 
    trans( "5pucksN.dat", "1824_N.dat", 0.39, 3.64, fitFile="1824.fitFile" ) 
    trans( "5pucksN.dat", "1825_N.dat", 0.39, 4.45, fitFile="1825.fitFile" ) 
    binit( "1818_N" )
    binit( "1822_N" )
    binit( "1823_N" )
    binit( "1824_N" )
    binit( "1825_N" )
    unpack( "1puckB.raw", "1puckB.dat" )
    trans( "1puckB.dat", "rex_B.dat", 0.39, 2.03, fitFile="rex.fitFile" ) 
    unpack( "5pucksO.raw", "5pucksO.dat" )
    trans( "5pucksO.dat", "1821_O.dat", 0.39, 1.22 ) 
    trans( "5pucksO.dat", "1820_O.dat", 0.39, 2.03 ) 
    trans( "5pucksO.dat", "1819_O.dat", 0.39, 2.82 ) 
    trans( "5pucksO.dat", "1827_O.dat", 0.39, 3.64 ) 
    trans( "5pucksO.dat", "1828_O.dat", 0.39, 4.45 ) 
    binit( "1821_O" )
    binit( "1820_O" )
    binit( "1819_O" )
    binit( "1827_O" )
    binit( "1828_O" )
    unpack( "5pucksA.raw", "5pucksA.dat" )
    trans( "5pucksA.dat", "1829_A.dat", 0.39, 1.22 ) 
    trans( "5pucksA.dat", "1830_A.dat", 0.39, 2.03 ) 
    trans( "5pucksA.dat", "1831_A.dat", 0.39, 2.82 ) 
    trans( "5pucksA.dat", "1832_A.dat", 0.39, 3.64 ) 
    trans( "5pucksA.dat", "1833_A.dat", 0.39, 4.45 ) 
    binit( "1829_A" )
    binit( "1830_A" )
    binit( "1831_A" )
    binit( "1832_A" )
    binit( "1833_A" )
    unpack( "5pucksB.raw", "5pucksB.dat" )
    trans( "5pucksB.dat", "1834_B.dat", 0.39, 1.22 ) 
    trans( "5pucksB.dat", "1835_B.dat", 0.39, 2.03 ) 
    trans( "5pucksB.dat", "1836_B.dat", 0.39, 2.82 ) 
    trans( "5pucksB.dat", "1837_B.dat", 0.39, 3.64 ) 
    trans( "5pucksB.dat", "1838_B.dat", 0.39, 4.45 ) 
    binit( "1834_B" )
    binit( "1835_B" )
    binit( "1836_B" )
    binit( "1837_B" )
    binit( "1838_B" )
def doit():
    unpack( "5pucksC.raw", "5pucksC.dat" )
    trans( "5pucksC.dat", "1834_C.dat", 0.39, 1.22 ) 
    trans( "5pucksC.dat", "1837_C.dat", 0.39, 2.03 ) 
    trans( "5pucksC.dat", "1838_C.dat", 0.39, 2.82 ) 
    trans( "5pucksC.dat", "1835_C.dat", 0.39, 3.64 ) 
    trans( "5pucksC.dat", "1836_C.dat", 0.39, 4.45 ) 
    binit( "1834_C" )
    binit( "1835_C" )
    binit( "1836_C" )
    binit( "1837_C" )
    binit( "1838_C" )
