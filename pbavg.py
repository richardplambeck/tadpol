# pbavg.py
# specialized routine to compute average 1mm mixer passbands
# note: all data files must have exactly the same number of points and freq range

import numpy

fin1 = open('/o/plambeck/rcvr/mxrtests/pblist.dat', "r")
flist = []
plist = []
tlist = []
nmxr = 0
npts = 0
for line1 in fin1 :
  if not line1.startswith("#") :
    a1 = line1.split()
    if len(a1) > 0 :
      fin2 = open( a1[0], "r" )
      nmxr = nmxr + 1
      npts = 0
      for line2 in fin2 :
        if not line2.startswith("#") :
          a2 = line2.split()
          flist.append( float(a2[0]) )
          plist.append( float(a2[1]) ) 
          tlist.append( float(a2[3]) ) 
          npts = npts + 1
      print "%s   npts = %d  t = %.2f" % (a1[0],npts, tlist[len(tlist)-100])
f = numpy.reshape(flist, (nmxr, -1) )
p = numpy.reshape(plist, (nmxr, -1) )
t = numpy.reshape(tlist, (nmxr, -1) )
favg = numpy.average(f, axis=0) 
pavg = numpy.average(p, axis=0) 
tavg = numpy.average(t, axis=0) 
pstd = numpy.std(p, axis=0)
tstd = numpy.std(t, axis=0)
fin1.close()
fin2.close()

fout = open("pbavg.dat", "w")
for fx,px,pstdx,tx,tstdx in zip (favg,pavg,pstd,tavg,tstd) :
  if (fx > 1.) :
    fout.write("%8.4f %8.4f %8.4f %8.4f %8.4f\n" % (fx,px,pstdx,tx,tstdx))
fout.close()
