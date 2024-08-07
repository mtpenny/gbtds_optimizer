import numpy as np
import pandas as pd
from pypolyclip import clip_multi
import matplotlib.pyplot as plt

"""
Convert a map into polygons
"""

filename = 'H22.6_eventrates.rates'
lbmap = pd.read_csv(filename,header=None,sep='\s+',names=['l','b','yield'])

#Figure out the grid -

#sort the map before converting it
lbmap.sort_values(['b','l'],inplace=True,ignore_index=True)

llist = lbmap.iloc[:,0]
blist = lbmap.iloc[:,1]
yieldlist = lbmap.iloc[:,2]

#Assumes a complete, evenly spaced map
lpix = np.unique(llist)
lspacing = np.average(lpix[1:]-lpix[:-1])
lorigin = lpix[0]
bpix = np.unique(blist)
bspacing = np.average(bpix[1:]-bpix[:-1])
borigin = bpix[0]
map_naxis = (bpix.shape[0],lpix.shape[0])

lmap = llist.to_numpy().reshape(map_naxis) 
bmap = blist.to_numpy().reshape(map_naxis)
yieldmap = yieldlist.to_numpy().reshape(map_naxis)

#lmesh,bmesh = = np.meshgrid(np.append(bpix,bpix[-1]+bspacing)-0.5*bspacing,
#                            np.append(lpix,lpix[-1]+lspacing)-0.5*lspacing,
#                            indexing='xy') 

lspike=38
bspike=33
yieldmap[bspike,lspike]=50000
print(lmap[bspike,lspike],bmap[bspike,lspike])
#yieldmap[map_naxis[0]-1,1]=20000
#print(lmap[map_naxis[0]-1,1],bmap[map_naxis[0]-1,1])
#yieldmap[1,map_naxis[1]-1]=30000
#print(lmap[1,map_naxis[1]-1],bmap[1,map_naxis[1]-1])
#yieldmap[map_naxis[0]-1,map_naxis[1]-1]=40000
#print(lmap[map_naxis[0]-1,map_naxis[1]-1],bmap[map_naxis[0]-1,map_naxis[1]-1])

#Load the field vertices
fieldsfile = open('layout_7f_3.chips','r')
chips = {}
cl = []
cb = []
cn = []

nv=0
nc=0
for line in fieldsfile:
    #print(line)
    if nv==0:
        vl = []
        vb = []
        vn = []
    
    if not line in ('\n', '\r\n'):
        lsplit = line.split(" ")
        vn.append(lsplit[0])
        vl.append(lsplit[5])
        vb.append(lsplit[6])
        nv+=1
    else:
        if nv>2:
            #Calculate the vertex coordinates on the pixel grid which has the
            #origin at the bottom left of the pixel grid
            cl.append(((np.array(vl).astype(float)-lorigin)/lspacing+0.5).tolist())
            cb.append(((np.array(vb).astype(float)-borigin)/bspacing+0.5).tolist())
            cn.append(vn[0])
            nv=0
            nc+=1
            #There are vertices, make a polygon

print(cl)
print(cb)
print(cn)
print(nc)

#chips_l = (np.array(cl)-lorigin)/lspacing
#chips_b = (np.array(cb)-borigin)/bspacing

lidx, bidx, area, slices = clip_multi(cl, cb, map_naxis)

# slices is a list of slice objects to link between the input polygons
# and the clipped pixel grid.
# clipped_l,clipped_b are the grid indices with overlapping pixels.
# area is the overlapping area on a given pixel.
# Each slice belongs to a polygon, and it has the number of elements corresponding
# to the number of pixels it covers
# Units of area are pixels
# The total yield is the sum of area * the yield/pixel of the pixel over slices



# the slices object can be used to get the area of each polygon
totalYield=0
totalArea=0
for i, s in enumerate(slices):
    print(lidx[s],lmap[bidx[s],lidx[s]])
    print(bidx[s],bmap[bidx[s],lidx[s]])
    print(area[s],yieldmap[bidx[s],lidx[s]])
    print(f'total area for polygon {i}={np.sum(area[s])}')
    print(f'total yield for polygon {i}={np.sum(area[s]*yieldmap[bidx[s],lidx[s]])}')
    totalArea += np.sum(area[s])
    totalYield += np.sum(area[s]*yieldmap[bidx[s],lidx[s]]*lspacing*bspacing)
print(totalArea,totalArea*lspacing*bspacing,totalYield)

plt.figure()
#rm = np.array(ratemap).reshape((lpix.shape[0],bpix.shape[0]))
plt.pcolormesh(lmap,bmap,yieldmap,shading='nearest')

for i,l in enumerate(cl):
    #plt.plot(lorigin+np.array(l)*lspacing,borigin+np.array(cb[i])*bspacing,'k-')
    plt.plot((np.array(l)-0.5)*lspacing+lorigin,
             (np.array(cb[i])-0.5)*bspacing+borigin,'k-')
#plt.gca().invert_xaxis()
plt.xlim([2.5,-2.5])
plt.ylim([-2.5,2.5])
plt.gca().set_aspect('equal')
plt.show()



#Test map and clipping with a map that is all zero except for small part(s)
