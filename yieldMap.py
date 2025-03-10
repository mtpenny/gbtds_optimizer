import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy


class yieldMap:

    def __init__(self,filename,read_csv_kwargs={},
                 spike=[],units='per_tile'):

        '''
        Constructor that takes a filename and pandas read_csv 
        keyword arguments to build a yieldMap object. The file 
        should have at least three columns: l and b coordinates 
        and a yield, and these columns should be named "l", "b",
        and "yield". The default arguments assume a 3-column file
        with columns, l, b, and yield with no header line. The
        map should be a regular grid with all values defined. It
        can be modified with "spike" - a list of 3-tuples with
        (l,b,yield) - which is useful for testing.
        '''

        self.filename = filename
        acceptable=False

        #Try to autodetect the separator if no sep keyword arg passed
        if not 'sep' in read_csv_kwargs:
            read_csv_kwargs['sep']=None
            read_csv_kwargs['engine']='python'

        #Read the map and account for different formats
        try:
            self.lbmap_orig = pd.read_csv(self.filename,
                                          on_bad_lines='error',
                                          dtype=float,**read_csv_kwargs)
            #usecols=['l','b','yield','alphaTexp','alphaCadence'],

            if {'l','b','yield','alphaTexp','alphaCadence'}.issubset(self.lbmap_orig.columns):
                acceptable=True
        except:
            pass


        if not acceptable:
            try:
                self.lbmap_orig = pd.read_csv(self.filename,
                                              usecols=['l','b','yield'],
                                              on_bad_lines='error',
                                              dtype=float,**read_csv_kwargs)

                if {'l','b','yield'}.issubset(self.lbmap_orig.columns):
                    acceptable=True                    
            except:
                pass

        if not acceptable:
            try:
                self.lbmap_orig = pd.read_csv(self.filename,header=None,
                                              names=['l','b','yield','alphaTexp',
                                                     'alphaCadence'],
                                              dtype=float,
                                              on_bad_lines='error',
                                              **read_csv_kwargs)

                if {'l','b','yield','alphaTexp','alphaCadence'}.issubset(self.lbmap_orig.columns):
                    acceptable=True
            except:
                pass


        if not acceptable:
            try:
                self.lbmap_orig = pd.read_csv(self.filename,header=None,
                                              names=['l','b','yield'],
                                              dtype=float,
                                              on_bad_lines='error',
                                              **read_csv_kwargs)
                
                if {'l','b','yield'}.issubset(self.lbmap_orig.columns):
                    acceptable=True                    
            except:
                raise RuntimeError('Error reading yieldMap file (%s). Try checking its formatting. It must have l,b, and yield columns, and can optionally have alphaTexp and alphaCadence columns too.' % (self.filename))

        #print(self.lbmap_orig)
        #print(self.lbmap_orig.shape)

        if units=='per_deg2' or units=='per-deg2' or units=='per_deg' or units=='per-deg':
            self.lpix = np.unique(self.lbmap_orig['l']) #Unique returns sorted list
            self.bpix = np.unique(self.lbmap_orig['b'])
            self.lspacing = np.average(self.lpix[1:]-self.lpix[:-1])
            self.bspacing = np.average(self.bpix[1:]-self.bpix[:-1])
            print("Assuming yield units are per deg2. Tile area is %g deg^2" % (self.lspacing*self.bspacing))
            self.lbmap_orig['yield'] *= self.bspacing * self.lspacing
        else:
            print("Assuming yield units are per tile")

        

        #Make a copy that will actually be used and modified
        self.lbmap_working = copy.deepcopy(self.lbmap_orig)


        self.processMap(spike)


    def processMap(self,spike=[]):

        '''
        Process a map from a pandas dataframe. The map is 
        converted into an integer grid and conversion methods 
        are supplied as part of the class. Spike is a list 
        of 3-tuples that can be used to modify the map using its
        integer coordinates - the main purpose is for testing.
        '''

        #sort the map before converting it
        self.lbmap = self.lbmap_working.sort_values(['b','l'],
                               ignore_index=True)

        #The map coordinates and values in 1d form
        
        self.llist = self.lbmap['l']
        self.blist = self.lbmap['b']
        self.yieldlist = self.lbmap['yield']


        #Assume complete, evenly spaced map & compute properties
        self.lpix = np.unique(self.llist)
        self.lspacing = np.average(self.lpix[1:]-self.lpix[:-1])
        self.lorigin = self.lpix[0]

        self.bpix = np.unique(self.blist)
        self.bspacing = np.average(self.bpix[1:]-self.bpix[:-1])
        self.borigin = self.bpix[0]

        self.map_naxis = (self.bpix.shape[0],self.lpix.shape[0])

        #The map coordiantes and values in 2d form
        self.lmap = self.llist.to_numpy().reshape(self.map_naxis)
        self.bmap = self.blist.to_numpy().reshape(self.map_naxis)
        self.yieldmap = self.yieldlist.to_numpy().reshape(self.map_naxis)
        self.covfac = np.zeros(self.yieldmap.shape)
        

        #Apply a "spike" to modify the map
        for sp in spike:
            print("Adding spike",sp)
            self.yieldmap[sp[1],sp[0]] = sp[2]


    def plotMap(self,ax=None,pcolormesh_kwargs={}):
        if ax is None:
            ret = plt.pcolormesh(self.lmap,self.bmap,self.yieldmap,
                                 shading='nearest',**pcolormesh_kwargs)
        else:
            ret = ax.pcolormesh(self.lmap,self.bmap,self.yieldmap,
                                shading='nearest',**pcolormesh_kwargs)
        return ret


    def l2x(self,l_vertices):
        '''
        Convert 2-d numpy array of l values to 2-d list of 
        map x pixel values
        '''
        return ((l_vertices.astype(float)-self.lorigin)/self.lspacing+0.5).tolist()


    def b2y(self,b_vertices):
        '''
        Convert 2-d numpy array of l values to 2-d list of 
        map x pixel values
        '''
        return ((b_vertices.astype(float)-self.borigin)/self.bspacing+0.5).tolist()
 
    def x2l(self,x):
        return (np.array(x)-0.5)*self.lspacing+self.lorigin

    def y2b(self,y):
        return (np.array(y)-0.5)*self.bspacing+self.borigin



        

        

        
