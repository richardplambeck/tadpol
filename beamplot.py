import numpy
import math
import matplotlib
import pickle
#matplotlib.use('GTK')
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.mlab import griddata
# from scipy.interpolate import griddata

# newStyle=False for data from Toki's LabView setup, newStyle=True for Campbell setup (sep 2018)
def readData( infile, newStyle=True ) :
  x = []
  y = []
  z = []
  phs = []
  fin = open( infile, "r" )
  for line in fin :
    if not line.startswith('#') :
      a = line.split()
      x.append(float(a[0]))
      y.append(float(a[1]))
      z.append(float(a[2]))
      if newStyle :
        phs.append(float(a[3]))
  fin.close()
  if newStyle :
    return numpy.array(x), numpy.array(y), numpy.array(z), numpy.array(phs)
  else :
    return numpy.array(x)/100000., numpy.array(y)/100000., numpy.array(z), numpy.zeros(len(x))

def makePlot( infile, x0=0., y0=0., distIn=6.2, Angle=False ) :
  xa,ya,z,phs = readData( infile )
  theta = numpy.arctan( numpy.sqrt(numpy.power(xa-x0,2)+numpy.power(ya-y0,2))/distIn )
  # zc = z/numpy.cos(theta)
  zc = z 
  if Angle :       # convert to angle
    x = numpy.degrees( numpy.arctan((xa-x0)/distIn))
    y = numpy.degrees( numpy.arctan((ya-y0)/distIn))
  else :            # convert to inches
    x = (xa - x0)
    y = (ya - y0)
  delta = (y[1]-y[0])/5.
  xi = numpy.arange( x.min(),x.max(),delta )
  yi = numpy.arange( y.min(),y.max(),delta )

  pyplot.figure()
  zi = griddata( x,y,zc, xi,yi, interp='linear' )
  cint = zi.max()/10.
  clevs = numpy.arange( cint, zi.max(), cint )
  halfpwr = [zi.max()/2.]
  pyplot.contour(xi,yi,zi, linewidths=1, levels=halfpwr, colors='k')
  pyplot.contour(xi,yi,zi, linewidths=0.5, levels=clevs, colors='k')
  pyplot.pcolormesh(xi,yi,zi,cmap=pyplot.get_cmap("rainbow"))
  pyplot.axes().set_aspect("equal")
  pyplot.title( infile )
  pyplot.grid()
  cbar = pyplot.colorbar()
  cbar.set_label("mV", rotation=90)
  pyplot.show()

# plot amp and phase of beamplots created with John's software (which dumps pkl file)
def ampPhsPlot( infile, x0=0., y0=0., distIn=6.2, Interpolate=False, ExtendPhase=False, Angle=False ) :

    xin = []
    yin = []
    ampin = []
    phsin = []

    fin = open( infile+'.pkl', 'rb')
    data = pickle.load( fin)
    fin.close()
    for i in range(0,len(data['phases'])) :
      xin.append( data['positions'][i][0] )
      yin.append( data['positions'][i][1] )
      ampin.append( data['amplitudes'][i] )
      phsin.append( data['phases'][i] )

    x = numpy.array(xin) - x0
    y = numpy.array(yin) - y0
    amp = numpy.array(ampin)
    phs = numpy.array(phsin)

    try :
      desc = data['description']
    except :
      desc = ''

    theta = numpy.arctan( numpy.sqrt(numpy.power(x-x0,2)+numpy.power(y-y0,2))/distIn )

  # create x and y arrays with 0,0 = center position
    if Angle :         # convert to angle
      x = numpy.degrees( numpy.arctan(x/distIn))
      y = numpy.degrees( numpy.arctan(y/distIn))

  # extend phases through +/- 180 degrees if desired
    if ExtendPhase :
      fout = open("extendphase.dat", "w")
      rindex = numpy.argsort( numpy.square(x)+numpy.square(y), axis=None )  
      for i in range(0,len(rindex)) :
        j = rindex[i]
        if i == 0:
          firstphs = phs[j]
          lastphs = phs[j]
          phs[j] = 0.
        else :
          oldphs = phs[j]
          while (phs[j]-lastphs) > 180. :
            phs[j] = phs[j] - 360.
          while (phs[j]-lastphs) < -180. :
            phs[j] = phs[j] + 360.
          # print "lastphs = %8.2f, oldphs = %8.2f, newphs = %8.2f" % (lastphs,oldphs,phs[j])
          fout.write("%9.3f  %9.3f  %9.3f\n" % (numpy.square(x[j]) + numpy.square(y[j]), oldphs, phs[j] ))
          lastphs = phs[j]
          phs[j] = phs[j] - firstphs   
      fout.close()
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()

  # if Interpolate, create interpolated arrays for smoother plots
    #if Interpolate :
    #  delta = (yin[1]-yin[0])/5.
    #else :
    #  delta = (yin[1]-yin[0])
    delta=0.25
    xi = numpy.arange( x.min(),x.max(),delta )
    yi = numpy.arange( y.min(),y.max(),delta )
    ampi = griddata( x,y, amp, xi,yi, interp='linear' )    # this is the mlab version
    phsi = griddata( x,y, phs, xi,yi, interp='linear' )
    #ampi = griddata( (x,y), amp, (xi,yi), method='linear' )
    #phsi = griddata( (x,y), phs, (xi,yi), method='linear' )

  # create side by side amp and phase plots
    fig, ax = pyplot.subplots( nrows=1, ncols=2, figsize=(14,5) )
  
    ampint = ampi.max()/20.
    amplevs = numpy.arange( ampint, ampi.max(), ampint )
    halfpwr = [ampi.max()/2.]
    ax[0].contour(xi,yi,ampi, linewidths=1, levels=halfpwr, colors='k')
    ax[0].contour(xi,yi,ampi, linewidths=0.5, levels=amplevs, colors='k')
    #ampplt = ax[0].pcolormesh(xi,yi,ampi,cmap=pyplot.get_cmap("rainbow"))
    ampplt = ax[0].imshow(ampi, extent=(xmin,xmax,ymin,ymax), cmap=pyplot.get_cmap("rainbow"), \
       aspect='equal', origin='lower' )
    ax[0].grid()
    cbar = fig.colorbar(ampplt, ax=ax[0])
    cbar.set_label("amplitude", rotation=90)

    #phsplt = ax[1].pcolormesh(xi,yi,phsi,cmap=pyplot.get_cmap("rainbow"))
    phslevs = numpy.arange( -3600., 3600.,180.)
    ax[1].contour(xi,yi,phsi, linewidths=0.5, levels=phslevs, colors='k')
    phsplt = ax[1].imshow(phsi, extent=(xmin,xmax,ymin,ymax), cmap=pyplot.get_cmap("rainbow"), \
       aspect='equal', origin='lower' )
    ax[1].set_aspect("equal")
    ax[1].grid()
    cbar = fig.colorbar(phsplt, ax=ax[1])
    cbar.set_label("phase", rotation=90)

    pyplot.suptitle( infile + "   " + desc )
    pyplot.show()

def unpack( infile ) :
    fin = open( infile+'.pkl', 'rb')
    data = pickle.load( fin)
    fin.close()
    fout = open( infile, "w" )
    for i in range(0,len(data['phases'])) :
      x = data['positions'][i][0]
      y = data['positions'][i][1]
      amp = data['amplitudes'][i]
      phs = data['phases'][i]
      fout.write("%9.5f  %9.5f  %9.5f  %9.5f\n" % (x,y,amp,phs ))
    fout.close()
   
   

#makePlot("20180208/16_beamMap_heterodyne_090GHz.txt", x0=10000., y0=25000., distIn=None )
#makePlot("20180320/6-20180320-intsph-90-corrhorn.txt", x0=0., y0=0., distIn=None )
# something is wrong with angle calc - beam turns square
# makePlot("20180320/6-20180320-intsph-90-corrhorn.txt", x0=0., y0=0., distIn=6.2 )
# makePlot("20180320/6-20180320-intsph-90-corrhorn.txt", x0=0., y0=0., distIn=None )
# makePlot("20180320/8-bolo-90.txt", x0=0., y0=0., distIn=None )
# makePlot("20180320/9-bolo-90_3secpause.txt", x0=0., y0=0., distIn=None )
