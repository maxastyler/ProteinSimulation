#!/usr/bin/env python3
import sys
import numpy as np
from pylab import *
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.spatial import Delaunay

# parameters
label='wt'
file='../../plumed/bulk/fes_137.7_avg.dat'
eunit=0.861718724 # 1/(2.479/298*139.5) , 1 kj/mol in kt
# set xran,yran
xvalues = np.arange(0.0,1.001,0.005)
yvalues = np.arange(0.0,1.001,0.005)
xi,yi = np.meshgrid(xvalues,yvalues)

nall=299
nbeta=133
nalpha=15
ncore=46
nctail=35
nab=nbeta+nalpha
noth=nall-nbeta-nalpha
nothnoc=nall-nctail-nbeta-nalpha
data = np.loadtxt(file)
print(data)
data = data[data[:,2]!=data[:,2]+1] # numbers different from nan,+inf,-inf
x =  1.0*data[:,0]/nab	# fraction of secondary structure contacts
y =  1.0*data[:,1]/nothnoc	# fraction of core contacts
z =  data[:,2]*eunit		# energy
min=np.nanmin(z)
z=z-min
# triangulate data
tri = Delaunay( np.array(list(zip(x,y))) )
# interpolate data
interp = CloughTocher2DInterpolator(tri,z)
zi = interp(xi,yi)
for i in range(len(xvalues)):
  for j in range(len(yvalues)):
    if(zi[i,j]<0.0): zi[i,j]=0.000001

# do the acutal plotting
levels = np.arange(0,18,2)
CSF = contourf(xi,yi,zi,levels)
CS = contour(xi,yi,zi,levels,colors=('k'))
bar = colorbar(CSF)
bar.set_label('energy (kT)')
#clabel(CS)
xlabel('alpha and beta contacts (fraction)')
ylabel('core contacts (fraction)')

# save
savefig('fes_bulk_wt_reweight.pdf')