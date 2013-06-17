# this routine is meant to convert cansniffer output to ascii file with times, integers
# to take the data:
#   - ssh to appropriate antenna computer (e.g. ssh c7)
#   - dump particular telem address to a file with cansniffer:
#       (cansniffer mid=0x10a mjd=true >> dump7)
#   - then run convert on this file to get list of times, integers

def convert( infile, outfile ) :
  fin = open( infile, "r" )
  fout = open( outfile, "w" )
  time0 = 0.
  for line in fin :
    a = line.split(",")
    b = a[0].split()
    time = float(b[1]) 
    if time0 == 0. :
      time0 = time
    c = a[5].split()
    i1 = int( c[3]+c[4], 16 )
    i2 = int( c[5]+c[6], 16  )
    i3 = int( c[7]+c[8], 16 )
    i4 = int( c[9]+c[10], 16 )
    print 1.e7*(time-time0), i1,i2,i3,i4
    fout.write("%10.3f %6d %6d %6d %6d\n" % (1.e5*(time-time0),i1,i2,i3,i4) )
  fin.close()
  fout.close()
