"""
Code to find the optimal telescope slew path through a small* number of fields given estimates of slew and settle times. The specific implementation allows for different slew times in different directions (Roman's WFI short axis, long axis, and diagonally).

*The code computes all possible cyclic permutations of the path through the fields, so compute and memory will scale as (N-1)! for large numbers of fields >~10

Matthew Penny (penny1@lsu.edu)

v0.1: Working
"""

import sys

if len(sys.argv) != 3 and len(sys.argv) != 5:
    print("")
    print("Usage: %s <fields> <slew-times (short axis)> {<slew-times (diagonal)> <slew-times (long-axis)>" % (sys.argv[0]))
    print("")
    print("For a given set of fields, find the path through the fields that gives the shortest total slew time for the given slew and settle times. Either a single slew time vs distance can be given, or three - one each for a short-axis, diagonal, and long-axis slew, respectively.")
    print("")
    print("Fields file: Name l(deg) b(deg)")
    print("Slew file: Distance(deg) Time(s)")
    sys.exit()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d
from itertools import permutations
from fields import fov

#Options
plotSlew = False
debug = True
    
#Load the fields file
fieldsFile = sys.argv[1]
fields = pd.read_csv(fieldsFile,sep=r'\s+',comment='#')
fieldNames = fields.iloc[:,0]
l = fields.iloc[:,1]
b = fields.iloc[:,2]
nfields = fields.shape[0]
if debug==True:
    print(l,b)
coords = SkyCoord(l,b,frame='galactic',unit='deg')

fix=True

if debug==True:
    print(coords)

#Load the slew time files
shortSlewFile = sys.argv[2]
diagSlewFile = shortSlewFile
longSlewFile = shortSlewFile

if len(sys.argv)==5:
    diagSlewFile = sys.argv[3]
    longSlewFile = sys.argv[4]

shortSlew = pd.read_csv(shortSlewFile,sep=r'\s+',header=None,comment='#')
shortSlewFn = interp1d(shortSlew.iloc[:,0],shortSlew.iloc[:,1],
                 fill_value="extrapolate")
diagSlew = pd.read_csv(diagSlewFile,sep=r'\s+',header=None,comment='#')
diagSlewFn = interp1d(diagSlew.iloc[:,0],diagSlew.iloc[:,1],
                fill_value="extrapolate")
longSlew = pd.read_csv(longSlewFile,sep=r'\s+',header=None,comment='#')
longSlewFn = interp1d(longSlew.iloc[:,0],longSlew.iloc[:,1],
                fill_value="extrapolate")

if plotSlew==True:
    plt.figure()
    x=np.arange(0.0,5.0,0.01)
    plt.plot(x,shortSlewFn(x),'k-')
    plt.plot(x,diagSlewFn(x),'r-')
    plt.plot(x,longSlewFn(x),'b-')
    plt.show()

#
# First task is to compute the matrix of all possible slew times between any pair of fields in
# the field pattern
#

    
#Construct the distance matrix of distances between two fields
sep = np.array([c.separation(coords) for c in coords])
angle = np.array([c.position_angle(coords).to(u.deg) for c in coords])

#Figure out if they are short, diagonal or long slews (0,1,2)
slewType = np.zeros(shape=sep.shape,dtype=int)

shortmask = ((angle>80.0) & (angle<100.0)) | ((angle>260.0) & (angle<280.0))
diagmask = ((angle>=10.0) & (angle<=80.0)) | ((angle>=100.0) & (angle<=170.0)) | ((angle>=190.0) & (angle<=260.0)) | ((angle>=280.0) & (angle<=350.0))
longmask = ((angle>-10.0) & (angle<10.0)) | ((angle>170.0) & (angle<190.0)) | ((angle>350) & (angle<370))

masks = [shortmask,diagmask,longmask]
slewType[diagmask] = 1
slewType[longmask] = 2
if debug==True:
    print(sep)
    print(angle)
    print(slewType)

#Compute the slew time
slewTimes = np.zeros(shape=sep.shape)
slewTimes[shortmask] = shortSlewFn(sep[shortmask])
slewTimes[diagmask] = diagSlewFn(sep[diagmask])
slewTimes[longmask] = longSlewFn(sep[longmask])

if debug==True:
    print(slewTimes)

#
# Next task is to compute all possible paths through each field once, then pick the best one
#

#Construct the permutations of the path through the fields 

def faster_cyc_permutations(m):
    """
    Permutations code copied from Daniel Giger's stack overflow reply
    https://stackoverflow.com/questions/64291076/generating-all-permutations-efficiently
    However, as we are only interested in cyclic permutations, we can cut down on the number 
    significantly by always starting/ending in the same place. This reduces the number to (m-1)!
    instead of m!, a saving of a factor of m in compute!
    """
    # empty() is fast because it does not initialize the values of the array
    # order='F' uses Fortran ordering, which makes accessing elements in the same column fast
    n=m-1
    perms = np.empty((np.math.factorial(n), n), dtype=np.uint8, order='F')
    perms[0, 0] = 0

    rows_to_copy = 1
    for i in range(1, n):
        perms[:rows_to_copy, i] = i
        for j in range(1, i + 1):
            start_row = rows_to_copy * j
            end_row = rows_to_copy * (j + 1)
            splitter = i - j
            perms[start_row: end_row, splitter] = i
            perms[start_row: end_row, :splitter] = perms[:rows_to_copy, :splitter]  # left side
            perms[start_row: end_row, splitter + 1:i + 1] = perms[:rows_to_copy, splitter:i]  # right side

        rows_to_copy *= i + 1

    return np.hstack((perms,n*np.ones(shape=(perms.shape[0],1),dtype=int)))


idxpaths = faster_cyc_permutations(nfields)
if debug==True:
    print(idxpaths.shape)
    print(idxpaths)
                    

if debug==True:
    np.set_printoptions(threshold=sys.maxsize)
    print(idxpaths)
    np.set_printoptions(threshold=20)


#Find the path that gives the shortest slewtimes 
minslew = 1.0e50
bestpath = np.array([])
for path in idxpaths:
    #roll shifts elements in array with wrapping
    #This also includes the return to start slew
    pathtime = np.sum(slewTimes[path,np.roll(path,1)])
    if fix:
        print(path,pathtime,np.roll(slewTimes[path,np.roll(path,1)],-1))
    
    if pathtime<minslew:
        minslew = pathtime
        bestpath = path

        
#Compute tables of cadence vs exposure time for choices of each
sdp = lambda x: "%.1f" % x
tdp = lambda x: "%.2f" % x

cadence1 = np.arange(10,21)
texp1 = (cadence1*60 - minslew)/nfields
alpha = 0.399 #Delta chi^2 power law scaling

#Yield relative to Penny+19 assuming new field(s) has no microlensing events
relyield1no = (1.0/(cadence1/15.0) * (texp1/46.8))**alpha
#Yield relative to Penny+19 assuming new field(s) has average microlensing planet detection rate of the original 7 fields
relyield1avg = relyield1no * (nfields/7.0)
#print(relyield1no)
#print(relyield1avg)

table1 = [["Choose cadence","",""],["Cadence(min)","Texp(s)", r"Rel. $1 M_{\oplus}$ yield"]]
for i,c in enumerate(cadence1):
    table1.append([sdp(c),sdp(texp1[i]),tdp(relyield1no[i]) + "-" + tdp(relyield1avg[i])])

#plt.figure()


texp2 = np.arange(14,25)*3.04
cadence2 = (texp2*nfields + minslew)/60.0

#Yield relative to Penny+19 assuming new field(s) has no microlensing events
relyield2no = (1.0/(cadence2/15.0) * (texp2/46.8))**alpha
#Yield relative to Penny+19 assuming new field(s) has average microlensing planet detection rate of the original 7 fields
relyield2avg = relyield2no * (nfields/7.0)
#print(relyield1no)
#print(relyield1avg)

table1 = [["Choose cadence","",""],["Cadence(min)","Texp(s)", r"Rel. $1 M_{\oplus}$ yield"]]
for i,c in enumerate(cadence1):
    table1.append([sdp(c),sdp(texp1[i]),tdp(relyield1no[i]) + "-" + tdp(relyield1avg[i])])


table2 = [["Choose Texp","",""],["Texp(s)","Cadence(min)", r"Rel. $1 M_{\oplus}$ yield"]]
for i,t in enumerate(texp2):
    table2.append([sdp(t),sdp(cadence2[i]),tdp(relyield2no[i]) + "-" + tdp(relyield2avg[i])])

print("Best combined slew time:",minslew)
print("Best average slew time per field:",minslew/nfields)
print("Best path:",bestpath)

#For plotting, need to add the first field to the
cycle = np.append(bestpath,bestpath[0])

#Plot the field layout and best path, and the tables
fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(10,6))

scas = fov('sca_layout.txt')


#Figure with path
ax1.plot(l,b,'ko',ms=5)
ax1.plot(l[cycle],b[cycle],'k-')
for i in range(len(l)): 
    ax1.annotate(fieldNames[i], (l[i], b[i] + 0.02))
    for chip in range(scas.nChips):
        ax1.plot(scas.delta_l[chip]+l[i],scas.delta_b[chip]+b[i],'k-',lw=0.5)

ax1.axis('square')
ax1.invert_xaxis()
ax1.set_aspect('equal')
ax1.set_xlabel("l (deg)")
ax1.set_ylabel("b (deg)")

#First table
ax2.axis("off")
ax2.axis("tight")
ax2.table(cellText=table1,loc="center",fontsize=24,colLoc='center')

#Second table
ax3.axis("off")
ax3.axis("tight")
ax3.table(cellText=table2,loc="center",fontsize=24,colLoc='center')

#Add a title
fig.text(0.5,0.9,fieldsFile + " - " + shortSlewFile,fontsize=14,ha='center')
fig.text(0.85,0.85,"Total slew time: %.1f s" % minslew,fontsize=12,ha='right')
fig.text(0.85,0.8,"Average slew time per field: %.1f s" % (minslew/nfields),fontsize=12,ha='right')
fig.text(0.85,0.75,"Best path: %s" % (fieldNames[bestpath].astype(str).str.cat(sep=' ')),fontsize=12,ha='right')
fig.tight_layout()

plt.savefig(fieldsFile + ".pdf")



                            
                   


