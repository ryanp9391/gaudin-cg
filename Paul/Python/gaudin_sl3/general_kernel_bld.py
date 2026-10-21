import sys
import general_kernel as gk
from fractions import Fraction as Fr
lam=tuple(int(c) for c in sys.argv[1])
tag=sys.argv[2]
th=gk.TH0
if len(sys.argv)>3:
    a,b=sys.argv[3].split(',')
    th=(Fr(a),Fr(b))
gk.get_point2(lam,tag,gk.ZETAS2[tag],th)
