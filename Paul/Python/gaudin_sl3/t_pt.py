import pickle, mpmath as mp
from fractions import Fraction as Fr
import general_kernel as gk
pt=pickle.load(open('general_kernel_pt_100_P_2.7_-3.5_80.pkl','rb'))
print(pt.keys()); print(pt['zeta'],pt['theta'],pt['d1'],pt['hw'])
s=pt['states'][0]; print(s.keys()); print(s['sector'],s['idx'],len(s['v']),len(s['dual']),s['degs'],s['ok'])
import glob
print(sorted(glob.glob('general_kernel_pt_*_P_*80.pkl')).__len__())
