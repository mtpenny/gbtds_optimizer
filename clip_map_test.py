import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fields import fov,fovHandler
from yieldMap import yieldMap



#Create a map
lspike=45
bspike=33
yspike=50000.0

filename = 'H22.6_eventrates.rates'
lbmap = yieldMap(filename)


#Load the field vertices

scafile = 'sca_layout.txt'
romanFoV = fov(scafile,unit='deg')

fieldsFile = 'field_layouts/layout_7f_3.centers'
fields = pd.read_csv(fieldsFile,sep='\s+',header=None,names=['field','l','b'])

#Create an fovHandler object and initialize it with the yield map, field
#of view and fields
handler = fovHandler()
handler.fromCentersChips(fields,romanFoV,lbmap)


plt.figure()
handler.plotFields()
lbmap.plotMap()
plt.xlim([2.5,-2.5])
plt.ylim([-2.5,2.5])
plt.gca().set_aspect('equal')




totalYield, totalAreaPix, totalArea = handler.computeYield()

print(totalYield,totalAreaPix,totalArea)



plt.show()
#Test map and clipping with a map that is all zero except for small part(s)


for offl in np.arange(-0.5,0.5,0.1):
    fieldsNew = fields.copy(deep=True)
    fieldsNew['l'] += offl
    handler.fromCentersChips(fieldsNew,romanFoV,lbmap)
    #lidx, bidx, area, slices = clip_multi(handler.lpix, handler.bpix, lbmap.map_naxis)
    #totalYield=0
    #for i, s in enumerate(slices):
    #    totalYield += np.sum(area[s]*lbmap.yieldmap[bidx[s],lidx[s]]*lbmap.lspacing*lbmap.bspacing)
    totalYield, totalAreaPix, totalArea = handler.computeYield()
    print(offl,totalYield)
