# print logarithmic contour intervals 

import sys
import getopt
import numpy

def main(argv) :
  a = numpy.zeros( 15 )
  print len(sys.argv)
  if len(sys.argv) < 3 :
    print "Usage: python logCont,py min factor"
    print "prints out 15 contours, min, factor*min, factor*factor*min, ..."    
  else :
    a[0] = float( sys.argv[1] )
    mult = float( sys.argv[2] )
    for n in range(1,15) :
      a[n] = a[n-1]*mult 
   
    print "levels %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f" \
        % (a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11])
   
if __name__ == "__main__" :
  main(sys.argv[1:])

