# readout.py
# these are routines related to Farzad's high frequency measurments

import numpy
import math
import sys
import time
import string
import pickle
import matplotlib
import datetime
matplotlib.use('GTKAgg')
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import pol
import scipy
from scipy.interpolate import splrep, splev

# lambda is alternate way of defining simple function
sin = lambda arg: numpy.sin(arg)
cos = lambda arg: numpy.cos(arg)
inv = lambda Mat: numpy.linalg.inv(Mat)
dr = lambda theta: math.pi*theta/180.
pow = lambda arg, index: numpy.power(arg, index)
sq = lambda arg: pow(arg, 2)


clight = 29.9792458 # speed of light, cm/nanosec

fMHz = numpy.arange(23.,29.01,0.01)
L = 4.7e-6
C = 8.2e-12
R = 10.
Zshunt = 1j*2.*math.pi*fMHz*1.e6*L + 1./(1j*2.*math.pi*fMHz*1.e6*C) + R
  # note that Zshunt is an array

# compute S11 and S21 for shunt resonator
# L is in Henrys, C in Farads, R in ohms
#
def shuntRes( fMHz, L, C, R ) :
    Zshunt = 1j*2.*math.pi*fMHz*1.e6*L + 1./(1j*2.*math.pi*fMHz*1.e6*C) + R
      # note that Zshunt is an array
    A = numpy.ones( len(fMHz), dtype=complex )
    B = numpy.zeros( len(fMHz), dtype=complex ) 
    C = 1./Zshunt
    D = numpy.ones( len(fMHz), dtype=complex )
    Z0 = 50.
    denom = A + B/Z0 + C*Z0 + D
    S11 = (A + B/Z0 - C*Z0 - D)/denom
    S21 = 2*(A*D - B*C)/denom
    S11log = 20.*numpy.log10(numpy.absolute(S11))
    S21log = 20.*numpy.log10(numpy.absolute(S21))
    return S11log, S21log

# reads data from network analyzer, returns fMHz,S arrays
def readData( infile ) :
    dataLine = False
    fMHz = []
    S = []
    fin = open( infile, "r" )
    for line in fin:
      if line.startswith("END") :
        dataLine = False
        fin.close()
        return numpy.array(fMHz), numpy.array(S)
      if dataLine :
        a = line.split(",")
        fMHz.append( float(a[0])/1.e6 )
        S.append( float(a[1]) )
      if line.startswith("BEGIN") :
        dataLine = True
      
      if line.startswith("BEGIN") :
        dataLines = True

def plotBoth( dataFile, L=4.7e-6, C=8.2e-12, R=10. ) :
    fMHz,S = readData( dataFile )
    S11,S21 = shuntRes( fMHz, L, C, R )
    pyplot.figure(0)
    pyplot.plot( fMHz, S21, color='red' )
    pyplot.plot( fMHz, S, color='blue' )
    pyplot.show()

Csrch = numpy.arange(8.2e-12,8.4e-12,.01e-12)
Lsrch = numpy.arange(4.4e-6,4.6e-6,.01e-6)
Rsrch = numpy.arange(9.,13,.05)
Osrch = numpy.arange(-.1,.11,.1)

def fit( dataFile ) :
   fMHz,S = readData(dataFile)
   vmin = 1.e6
   for C in Csrch :
     for L in Lsrch :
       for R in Rsrch :
         for O in Osrch :
           S11fit,S21fit = shuntRes( fMHz, L=L, C=C, R=R )
           v = numpy.var( S21fit-S+O )
           if v < vmin :
             vmin = v
             Lbest = L
             Cbest = C
             Rbest = R
             Obest = O
   print Lbest, Cbest, Rbest, Obest
   plotBoth( dataFile, L=Lbest, C=Cbest, R=Rbest )           
