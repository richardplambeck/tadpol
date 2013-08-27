# makes alphabetical list of all my coauthors from a bib file
# this is supposed to ease to job of listing all collaborators in past 48 months
# to use: copy relevant bib items to /o/plambeck/CV/Last48.bib
#    

infile = "/o/plambeck/CV/Plambeck_Publications.bib"
infile = "/o/plambeck/CV/Last48.bib"
fin = open( infile, "r")
authlist = False
alist = "" 
fulllist = ""

for line in fin:
  a = line.split()
  if len(a) > 0 :
    if "=" in line :
      if a[0] == "author" :
        authlist = True
        begin = line.find('{') + 1
        tmp = line[begin:].translate( None, '{}\n\t' ).rstrip( '},' ).replace( "and", ";" ) 
        alist = tmp
        
      else :
        authlist = False
        if len(alist) > 0 :
          print ""
          print alist
          if len(fulllist) > 0 :
            fulllist = fulllist + " ; " + alist
          else :
            fulllist = alist
        alist = ""

    else :
      if authlist :
        tmp = line.translate( None, '{}\n\t' ).rstrip( '},' ).replace( "and", ";" ) 
        alist = alist + tmp

fin.close()
print fulllist
print " "
auth = fulllist.split(";")
al = []
for au in auth :
  if not au.strip() in al :
    al.append(au.strip())
fout = open( "Collaborators.txt", "w" )
for x in sorted(al) :
  print x
  fout.write("%s\n" % x)
fout.close()

      

