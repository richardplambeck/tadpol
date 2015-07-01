import numpy
import math
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import os

# create trcvr[nLOfreq, nIFfreq] array

nIFfreq = 133

# splist is a list of lists, each containing file names created 
#   using bimaRx sptrace=1 > file
# first file MUST be room temp, 2nd MUST be LN2
 
def calcT( splist, plottitle=None, pdftitle='junk.pdf', dataOut='junk.dat', cutoff=0.08 ) :
  tamb = 295.
  tcold = 77.
  color = ['red','blue','teal','darkgreen','magenta','navy']
  f = []       # list of frequencies in GHz
  p = []       # list of powers in mW
  nf = []      # list containing number of files of each type
  fcolor = []  # list with color code for each file (red for amb, etc.)
  nfiles = 0   # total number of files of all types 

  # open files one by one, append freq and pwr to lists
  for nfl,filelist in enumerate(splist) :
    i = 0
    for file in filelist :
      if os.path.isfile("./"+file) :                # check if file exists
        nfiles = nfiles + 1
        fcolor.append( color[nfl] )
        i = i + 1
        fin = open( file, "r" )
        for line in fin :
          if not line.startswith("#") :
            a = line.split()
            f.append( float(a[0])/1000. )			# convert freq to GHz
            p.append( pow(10., float(a[1])/10.) )   # convert dBm to mW
        fin.close()
      else :
        print "warning: file ./%s does not exist" % file
    nf.append( i )

  # now convert powers to matrix; keep only 1st column of freq matrix
  fGHz = numpy.reshape( numpy.array(f), [nfiles,-1])[1]
  pwr =  numpy.reshape( numpy.array(p), [nfiles,-1])
  msk = numpy.zeros( (len(fGHz)), dtype=bool )
     # msk will be used to mask off outliers

  # make second pass through list to average like powers together
  pavg = []
  n1 = 0
  #pT = pyplot.subplot(1,1,1)
  for i,filelist in enumerate(splist) :
    n2 = n1 + nf[i]
    pavg.append( numpy.average( pwr[n1:n2], axis=0 ))
       # average of all the measurements at same freq
    pvar = numpy.std( pwr[n1:n2], axis=0 )
    pmean = numpy.average( pwr[n1:n2], axis=0 )
    varnorm = pvar/pmean 
    for j in range(0, len(varnorm)) :
      if fGHz[j] < 0.2 :
        msk[j] = True
      elif varnorm[j] > cutoff :
        print "... %.4f  %.3f  %s" % (fGHz[j],varnorm[j],filelist)
        if (i < 2 ) :
          msk[j] = True 
    #pT.plot( fGHz, varnorm )
    n1 = n2
  #pT.grid( True, linewidth=0.1, color="0.05" )  # color=0.05 is a light gray
  #pyplot.title( plottitle )
  # pyplot.show()     # plot normalized deviations

  # if p1 and p2 arrays exist, then compute Trcvr, gain, and plot all temps
  if len(splist) > 1 :
    gain = (pavg[0]-pavg[1])/(tamb-tcold)
    trcvr = pavg[1]/gain - tcold
    pyplot.ioff()
    pyplot.clf()
    pp = PdfPages( pdftitle )
    pT = pyplot.subplot(1,1,1)
    pT.axis( [0.2,10.,50.,330.] )
    # pT.plot( fGHz, trcvr )

    # plot average scans
    for i,pav in enumerate(pavg) :
      T = pav/gain - trcvr
      Tmasked = numpy.ma.masked_where( msk,  T )
      fmasked = numpy.ma.masked_where( msk, fGHz )
      # pT.plot( fGHz, T, linewidth=1 )   # this is the unmasked version
      Tavg = numpy.ma.mean( Tmasked )
      pT.plot( fmasked, Tmasked, linewidth=1, label="avg = %.2f K" % Tavg, color=color[i] )
      
    # now (if gain and Trcvr were calculated), plot every individual file as temp
    for i in range(0,nfiles) :
      T = pwr[i]/gain - trcvr
      pT.plot( fGHz, T, ',', color=fcolor[i] )
    pT.grid( True, linewidth=0.1, color="0.05" )  # color=0.05 is a light gray
    pT.legend( loc=1, prop={'size':6} )
    pyplot.title( plottitle )
    pyplot.savefig( pp, format='pdf', bbox_inches='tight' )
    pp.close()
    # pyplot.show() 

    # 7/1/15 - for debugging purposes, plot individual files also
    fout = open("junk.dat", "w")
    for i in range( 0, len(fGHz) ) :
      fout.write("%8.5f" % fGHz[i] )
      for j in range(0,nfiles) :
        fout.write("%8.2f" % (pwr[j][i]/gain[i] - trcvr[i]) )
      fout.write("\n")
    fout.close()

    # dump freq, T3, T4, T5, T6 to output file
    fout = open( dataOut, "w")
    fout.write("# %s\n" % plottitle )
    for i in range( 0, len(fGHz) ) :
      if not msk[i] :
        fout.write( "%8.5f  %8.3f  %8.3f  %8.3f  %8.3f\n" % \
           ( fGHz[i], pavg[2][i]/gain[i] - trcvr[i], pavg[3][i]/gain[i] - trcvr[i],
           pavg[4][i]/gain[i] - trcvr[i], pavg[5][i]/gain[i] - trcvr[i] ) )
    fout.close() 
   
      

def convertTrace( infile, outfile ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  for line in fin :
    if line.startswith("#") :
      fout.write( "%s\n" % line )
    else :
      a = line.split()
      fGHz = float(a[0])/1000.
      pwr = pow( 10., float(a[1])/10. )
      if fGHz > 1. :
        fout.write( "%10.4f  %8.5e\n" % (fGHz, pwr) )
  fin.close()
  fout.close()

def pbavg(start,stop,interval) :
  nLOfreq = (stop - start)/interval + 1
  print "nLOfreq =", nLOfreq
  trcvr = numpy.ones([nLOfreq,nIFfreq])
  freqIF = numpy.ones([nIFfreq])				       

  nLO = -1
  for suffix in range (start, stop+interval, interval) :
    nLO += 1
    fin = open("pb." + str(suffix), "r")
    nIF = -1					
    for line in fin:
      if not line.startswith("#") :
        a = line.split()
	nIF += 1
	freqIF[nIF] = float(a[0])	
        trcvr[nLO,nIF] = float(a[3])
    fin.close()

  tavg = numpy.ones([nIFfreq])
  rms = numpy.ones([nIFfreq])
  fout = open("pb.avg", "w")
  for nIF in range(0, nIFfreq) :
    sum = 0.
    for nLO in range(0, nLOfreq) :
      sum += trcvr[nLO,nIF]
    tavg[nIF] = sum/float(nLOfreq)
    sum = 0.
    for nLO in range(0, nLOfreq) :
      sum += (trcvr[nLO,nIF] - tavg[nIF])**2
    rms[nIF] = math.sqrt(sum/(nIFfreq - 1))
    print freqIF[nIF], tavg[nIF], rms[nIF]
    sum = 0.
    npts = 0	    
    for nLO in range(0, nLOfreq) :
      if (math.fabs((trcvr[nLO,nIF] - tavg[nIF])/rms[nIF]) > 2) :
	print nLO, freqIF[nIF], tavg[nIF], trcvr[nLO,nIF]
      else :
        sum += trcvr[nLO,nIF]
	npts += 1	 	     
    tavg[nIF] = sum/float(npts)
    fout.write("%10.5f %10.5f %8.5f %3d\n" % (freqIF[nIF], tavg[nIF], rms[nIF], npts))
  fout.close()
								      
  fout2 = open("pball","w")
  nLO = -1
  for suffix in range (start, stop+interval, interval) :
    nLO += 1
    fout = open("pb." + str(suffix) + ".norm", "w")
    for nIF in range (0, nIFfreq) :
      fout.write("%10.5f %10.4f\n" % (freqIF[nIF], trcvr[nLO,nIF]/tavg[nIF]) )
      fout2.write("%4d  %10.5f  %10.4f\n" % (nLO, freqIF[nIF], trcvr[nLO,nIF]/tavg[nIF]))
    fout.close()	  	
    fout2.write("\n")
  fout2.close()
