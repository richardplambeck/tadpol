# TbUranus.py
# return brightness temp of Uranus at frequency fGHz
# encodes model of Griffin & Orton 1993

import math

a0 = -795.694
a1 = 845.179
a2 = -288.946
a3 = 35.200

fGHz = 226.3
lambda_mm = 299.8/fGHz
print "lamba = %.3f mm" % (lambda_mm)
lda = 1000. * lambda_mm

T = a0 + a1 * math.log10(lda) + a2 * pow( math.log10(lda), 2 ) + a3 * pow( math.log10(lda), 3 ) 

print "Tb = %.2f" % (T)
