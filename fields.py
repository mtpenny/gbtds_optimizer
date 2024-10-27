import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pypolyclip import clip_multi
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d
from itertools import permutations
import sys
import copy

class fov:

    '''
    A class for holding a single field of view with vertices of
    chips/SCAs measured relative to the field of view center in
    units of degrees
    '''

    def __init__(self,filename,unit='deg',name_col=0,delta_l_col=1,
                 delta_b_col=2):

        try:
            f = open(filename,'r')
        except:
            raise RuntimeError('Error reading SCA file (%s)' % (filename))
        scale=1.0
        if unit=='deg' or unit=='degree' or unit=='degrees':
            scale=1.0
        if unit=='arcsec':
            scale=1.0/3600.0
        if unit=='arcmin':
            scale=1.0/60.0
        if unit=='rad' or unit=='radian':
            scale = 180.0/np.pi

        self.chip = []
        self.delta_l = []
        self.delta_b = []
        

        nv = 0
        self.nChips=0

        for line in f:
            if nv==0:
                vl = []
                vb = []
                vn = []
    
            if not line in ('\n', '\r\n'):
                lsplit = line.split(" ")
                vn.append(lsplit[name_col])
                vl.append(lsplit[delta_l_col])
                vb.append(lsplit[delta_b_col])
                nv+=1
            else:
                if nv>2:
                    #There are vertices, make a polygon

                    #Calculate the vertex coordinates on the
                    #pixel grid which has the origin at the
                    #bottom left of the pixel grid
                    
                    self.delta_l.append(np.array(vl).astype(float)*scale)
                    self.delta_b.append(np.array(vb).astype(float)*scale)
                    self.chip.append(vn[0])
                    nv=0
                    self.nChips+=1
                    

class fovHandler:

    '''
    Class for building/handling the vertices of chips or SCAs
    '''
    
    def __init__(self):
        '''
        Assumes vertices is a 4-tuple of 2-d numpy arrays
        '''
        #(self.field,self.chip,self.l,self.b) = vertices
        pass
        


    #@classmethod
    #def fromCentersFileChipsFile(self,centersFile,chipsFile):
    #    '''
    #    Build a set of vertices from the

    def fromCentersChips(self,centers,chips,yieldMap,
                         centers_kwargs={'sep':'\s+'},
                         debug=False):
        '''
        Build vertices from a centers object (representing field
        centers) and a chips object (representing the corners of
        the chips/SCAs in a field of view). Centers should be a
        pandas dataframe with columns 'field','l','b',{'fixed'} chips 
        should be a pandas data frame with columns 'chip','delta_l', 'delta_b'
        '''

        self.field = []
        self.chip = []
        self.lpix = []
        self.bpix = []
        self.yieldMap = yieldMap
        self.debug = debug

        if isinstance(centers,str):
            try:
                centers = pd.read_csv(centers,**centers_kwargs)
            except:
                raise RuntimeError('Error reading field centers file (%s)' % (centers))

        if not set(['l','b','field','fixed']).issubset(centers.columns):
            if set(['l','b','field']).issubset(centers.columns):
                centers['fixed']=0
            else:
                raise "centers is not a dataframe containing at least the columns 'l','b', and 'field'"

        #print(centers)

        for idx,row in centers.iterrows():
            for i in range(chips.nChips):
                self.lpix.append(yieldMap.l2x(row['l']
                                           +chips.delta_l[i]))
                self.bpix.append(yieldMap.b2y(row['b']
                                           +chips.delta_b[i]))
            self.chip.append(chips.chip[i])
            self.field.append(row['field'])

        self.lpix = np.array(self.lpix)
        self.bpix = np.array(self.bpix)


    def plotFields(self,ax=None,
                   plot_kwargs={'color':'k','linestyle':'-'}):
        for i,l in enumerate(self.lpix):
            if ax is None:
                #print(self.chip[i],self.lpix[i],self.bpix[i])
                plt.plot(self.yieldMap.x2l(self.lpix[i]),
                         self.yieldMap.y2b(self.bpix[i]),
                         **plot_kwargs)
            else:
                ax.plot(self.yieldMap.x2l(self.lpix[i]),
                        self.yieldMap.y2b(self.bpix[i]),
                        **plot_kwargs)


    def scaleMap(self,Cadence,C0,alphaC,Texp,Texp0,alphaTexp):
        if (np.isscalar(alphaC) and alphaC==0) and (np.isscalar(alphaTexp) and alphaTexp ==0):
            pass
        else:
            self.yieldMap.lbmap_working['yield'] = self.yieldMap.lbmap_orig['yield'] * \
                (Cadence/C0)**alphaC * (Texp/Texp0)**alphaTexp
            self.yieldMap.processMap()

    def computeYield(self):

        '''
        Compute the yield of the fovHandler's current layout and map.
        Returns totalYield, totalArea(map pixels), totalArea(deg**2)
        '''

        debug=False

        try:
            self.yieldMap
        except:
            raise NameError("Error: fovHandler.yieldMap Yield map not defined, likely because the fovHandler class has not been initialized.")

        try:
            self.lpix
        except:
            raise NameError("fovHander.lpix not defined, likely because the fovHandler class has not been initialized.")

        ym = self.yieldMap

        lidx, bidx, area, slices = clip_multi(self.lpix, self.bpix, ym.map_naxis)


        # slices is a list of slice objects to link between the input
        # polygons and the clipped pixel grid.
        # clipped_l,clipped_b are the grid indices with overlapping pixels.
        # area is the overlapping area on a given pixel.
        # Each slice belongs to a polygon, and it has the number of elements
        # corresponding to the number of pixels it covers
        # Units of area are pixels
        # The total yield is the sum of area * the yield/pixel of the
        # pixel over slices

        # the slices object can be used to get the area of each polygon
        totalYield=0
        totalArea=0

        if debug:
            print(lidx)
            print(bidx)

        for i, s in enumerate(slices):
            if debug:
                print(s)
                print(lidx[s].tolist())
                print(type(bidx[s]))
                print(lidx[s],ym.lmap[bidx[s],lidx[s]])
                print(bidx[s],ym.bmap[bidx[s],lidx[s]])
                print(area[s],ym.yieldmap[bidx[s],lidx[s]])
                print(f'total area for polygon {i}={np.sum(area[s])}')
                print(f'total yield for polygon {i}={np.sum(area[s]*ym.yieldmap[bidx[s],lidx[s]])}')
            totalArea += np.sum(area[s])
            totalYield += np.sum(
                area[s]*ym.yieldmap[bidx[s],lidx[s]])
        if self.debug:
            print(totalArea,totalArea*ym.lspacing*ym.bspacing,totalYield)

        return totalYield,totalArea,totalArea*ym.lspacing*ym.bspacing


class slewOptimizer:

    def __init__(self, shortSFile, diagSFile=None, longSFile=None, debug=False):

        #Load the slew time files
        self.shortSlewFile = shortSFile
        self.diagSlewFile = self.shortSlewFile
        self.longSlewFile = self.shortSlewFile
        self.debug=debug

        if diagSFile is not None:
            self.diagSlewFile = diagSFile
        if longSFile is not None:
            self.longSlewFile = longSFile

        self.shortSlew = pd.read_csv(self.shortSlewFile,sep='\s+',header=None,
                                     comment='#')
        self.shortSlewFn = interp1d(self.shortSlew.iloc[:,0],self.shortSlew.iloc[:,1],
                                    fill_value="extrapolate")
        self.diagSlew = pd.read_csv(self.diagSlewFile,sep='\s+',header=None,
                                    comment='#')
        self.diagSlewFn = interp1d(self.diagSlew.iloc[:,0],self.diagSlew.iloc[:,1],
                                   fill_value="extrapolate")
        self.longSlew = pd.read_csv(self.longSlewFile,sep='\s+',header=None,
                                    comment='#')
        self.longSlewFn = interp1d(self.longSlew.iloc[:,0],self.longSlew.iloc[:,1],
                              fill_value="extrapolate")

    def faster_cyc_permutations(self, m):
        """
        Permutations code copied from Daniel Giger's stack overflow reply
        https://stackoverflow.com/questions/64291076/generating-all-permutations-efficiently
        However, as we are only interested in cyclic permutations, we can cut down on the number 
        significantly by always starting/ending in the same place. This reduces the number to (m-1)!
        instead of m!, a saving of a factor of m in compute!
        """
        # empty() is fast because it does not initialize the values of the array
        # order='F' uses Fortran ordering, which makes accessing elements in the same
        # column fast
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

            
    def optimizePath(self,centers,fixPath=False):
        '''
        Find the optimal path through the field centers and return the total 
        overhead and the best path through the fields. if fixPath=True, just 
        compute the overhead for the path through the field in the order given.
        '''

        debug=False

        fieldNames = centers['field']
        l = centers['l']
        b = centers['b']
        #fixed = centers['fixed'].str.contains('fixed')
        fixed = centers['fixed']
        nfields = centers.shape[0]
        if self.debug==True:
            print('l:',' '.join(l.astype(str)))
            print('b:',' '.join(b.astype(str)))
        coords = SkyCoord(l,b,frame='galactic',unit='deg')

        #Construct the distance matrix of distances between two fields
        sep = np.array([c.separation(coords) for c in coords])
        angle = np.array([c.position_angle(coords).to(u.deg) for c in coords])

        #Figure out if they are short, diagonal or long slews (0,1,2)
        slewType = np.zeros(shape=sep.shape,dtype=int)

        shortmask = ((angle>80.0) & (angle<100.0)) | ((angle>260.0) & (angle<280.0))
        diagmask = ((angle>=10.0) & (angle<=80.0)) | \
        ((angle>=100.0) & (angle<=170.0)) | ((angle>=190.0) & (angle<=260.0)) | \
        ((angle>=280.0) & (angle<=350.0))
        longmask = ((angle>-10.0) & (angle<10.0)) | ((angle>170.0) & (angle<190.0)) | \
        ((angle>350) & (angle<370))

        masks = [shortmask,diagmask,longmask]
        slewType[diagmask] = 1
        slewType[longmask] = 2
        if debug==True:
            print(sep)
            print(angle)
            print(slewType)

        #Compute the slew time
        slewTimes = np.zeros(shape=sep.shape)
        slewTimes[shortmask] = self.shortSlewFn(sep[shortmask])
        slewTimes[diagmask] = self.diagSlewFn(sep[diagmask])
        slewTimes[longmask] = self.longSlewFn(sep[longmask])

        if debug==True:
            print(slewTimes)

        # Next task is to compute all possible paths through each field once, then
        # pick the best one


        #Construct the permutations of the path through the fields 

        if fixPath:
            idxpaths = np.array([range(nfields)])
        else:
            idxpaths = self.faster_cyc_permutations(nfields)

        if self.debug:
            print(idxpaths.shape)
            print(idxpaths)

        #if self.debug:
        #    np.set_printoptions(threshold=sys.maxsize)
        #    print(idxpaths)
        #    np.set_printoptions(threshold=20)


        #Find the path that gives the shortest slewtimes 
        minslew = 1.0e50
        bestpath = np.array([])
        for path in idxpaths:
            #roll shifts elements in array with wrapping
            #This also includes the return to start slew
            pathtime = np.sum(slewTimes[path,np.roll(path,1)])
            if self.debug:
                print(' '.join(fieldNames[path]),pathtime,np.roll(slewTimes[path,np.roll(path,1)],-1))
    
            if pathtime<minslew:
                minslew = pathtime
                bestpath = path

        return minslew,bestpath
        
        

                
        
