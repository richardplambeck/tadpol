# calcbreaks.py interprets schedule file

def doit( infile, outfile ) :
  f = open(infile, 'r')
  fout = open(outfile, 'w')
  reached_timedata = False

  for line in f:

	if (len(line.strip().split()) > 8):

        # if we've read at least one line, then time2 (previous stop time) exists
		if (reached_timedata == True):
		  time1 = line.strip().split()[5]
		  time1_bits = time1.split(":")
		  time1h = int(time1_bits[0])
		  time1m = int(time1_bits[1])
		  time1s = int(time1_bits[2])
		  timediff = (time1h - time2h) * 3600 + (time1m - time2m) * 60 + (time1s - time2s)
		  fout.write("gap %.1f min\n" % (timediff/60.) )

        # for first or any other line, save stop time
		time2 = line.strip().split()[8]
		time2_bits = time2.split(":")
		time2h = int(time2_bits[0])
		time2m = int(time2_bits[1])
		time2s = int(time2_bits[2])
		reached_timedata = True

		print line.strip()
        fout.write("%s\n" % line.strip() )
  f.close()
  fout.close()

