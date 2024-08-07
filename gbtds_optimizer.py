import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import pickle

from fields import fov,fovHandler,slewOptimizer
from yieldMap import yieldMap

#Create a map
#lspike=45
#bspike=33
#yspike=50000.0


#Plan of attack

#Fields can be fixed or mobile. All mobile fields will move together.
#At the moment, no test is done to handle overlapping fields, and yield will
#be double counted in these areas - this can be avoided by constraining the
#movement of fields to areas that will not cause overlap.

#The CCS definition committee will need to decide between many science cases
#so the most useful output is a contour map of the yield with a maximum
#indicated.

parser = argparse.ArgumentParser(description="Optimizer for the Roman Galactic Bulge Time Domain Survey")
parser.add_argument('YieldMap_Filename',
                    help='File name for the yield map')
parser.add_argument('Cadence0',type=float,
                    help='Cadence the map was computed for, in minutes')
parser.add_argument('texp0',type=float,
                    help='Exposure time the map was computed for, in seconds')
parser.add_argument('--Cadence_Bounds',nargs=2,default=[5.0,16.0],
                    help='Bounds on the cadence to be considered')
parser.add_argument('--NRead_Bounds',nargs=2,default=[6,30],
                    help='Bounds on the number of reads in an exposure')
#parser.add_argument('--dYdCadence',default=1.0,
#                    help='Map filename with the same shape as yieldMap or a numerical value of dY/dCadence')
#parser.add_argument('--dYdtexp',default=1.0,
#                    help='Map filename with the same shape as yieldMap or a numerical value of dY/dtexp')
parser.add_argument('--alphaCadence',default=None,type=str,
                    help='The power law slope of the yield as a function of Cadence, with input as a map filename with the same shape as yieldMap or a single numerical value')
parser.add_argument('--alphaTexp',default=None,type=str,
                    help='The power law slope of the yield as a function of Texp, with input as a map filename with the same shape as yieldMap or a single numerical value')
parser.add_argument('Fields_Filename',
                    help='File containing the field specifications (columns: name l b fixed?)')
parser.add_argument('--SCA_Filename',default='sca_layout.txt',
                    help='Filename of the sca vertices')
parser.add_argument('--roll',default=0.0,choices=[0.0,180.0],type=float,
                    help='Roll angle of the SCA layout. Options are 0 and 180')
parser.add_argument('--lrange',nargs=2,default=[3.0,-3.0],type=float,
                    help='Range of l in the grid')
parser.add_argument('--lstep',default=0.01,type=float,
                    help='Stepsize of the l grid')
parser.add_argument('--brange',nargs=2,default=[-3.0,0.0],type=float,
                    help='Range of b in the grid')
parser.add_argument('--bstep',default=0.01,type=float,
                    help='Stepsize of the b grid')
parser.add_argument('--readTime',default=3.04,type=float,
                    help='Time between up-the-ramp reads of a pixel in seconds')

parser.add_argument('--SlewRatesShortFoV_Filename',nargs=1,
                    default='slew_times_withResetReference_McEnery05232024.txt',
                    help='Filename for the short FoV slew and settle times')
parser.add_argument('--SlewRatesDiagonalFoV_Filename',nargs=1,default=None,
                    help='Filename for the diagonal FoV slew and settle times')
parser.add_argument('--SlewRatesLongFoV_Filename',nargs=1,default=None,
                    help='Filename for the long FoV slew and settle times')

parser.add_argument('--testyield',default=False,
                    help='Compute the yield for the input fields as a test')
parser.add_argument('--testplot',default=False,
                    help='Plot the layout of the input fields on the map a test')

parser.add_argument('--output_root',default='test',
                    help='Filename root for output')

args = parser.parse_args()

lstep = args.lstep; bstep = args.bstep;
lrange = args.lrange; brange = args.brange

#Load the yieldMap and its derivatives
ym = yieldMap(args.YieldMap_Filename)

#Load the power laws if needed
#Process some arguments
if args.alphaCadence is not None:
    try:
        alphaC = float(args.alphaCadence)
    except ValueError:
        try:
            alphaCadence = args.alphaCadence
            if args.alphaCadence == 'same':
                alphaCadence=args.YieldMap_Filename
            alphaC = pd.read_csv(alphaCadence,sep='\s+',usecols=['l','b','alphaCadence'])['alphaCadence']
        except:
            raise RuntimeError('Error reading alphaCadence file (%s)' % (alphaCadence))

if args.alphaTexp is not None:
    try:
        alphaT = float(args.alphaTexp)
    except ValueError:
        try:
            alphaTexp = args.alphaTexp
            if args.alphaTexp == 'same':
                alphaTexp=args.YieldMap_Filename
            alphaT = pd.read_csv(alphaTexp,sep='\s+',usecols=['l','b','alphaTexp'])['alphaTexp']
        except:
            raise RuntimeError('Error reading alphaTexp file (%s)' % (alphaTexp))



#Load the field vertices
romanFoV = fov(args.SCA_Filename,unit='deg')

try:
    fields = pd.read_csv(args.Fields_Filename,sep='\s+',header=None,names=['field','l','b','fixed'],
                         dtype={'field':str,'l':float,'b':float,'fixed':str})
except:
    try:
        fields = pd.read_csv(args.Fields_Filename,sep='\s+')
    except:
        raise RuntimeError('Error reading fields file (%s) - file does not exist or is in an incorrect format. It should containt the named columns ["field","l","b","fixed"] or 4 unnamed columns' % (args.Fields_Filename))

print("Using fields:")
print(fields)



allfixed = False
allfree = True
fixedMask = fields['fixed'].str.contains('fixed')
freeMask = np.logical_not(fixedMask)
nfixed = fixedMask.sum()
if nfixed==fields.shape[0]:
    allfixed = True
if nfixed>0:
    allfree = False

#Create an fovHandler object and initialize it with the yield map, field
#of view and fields
handler = fovHandler()
handler.fromCentersChips(fields,romanFoV,ym)

#Create a slewOptimizer object and initialize it with the slew times
slewopt = slewOptimizer(args.SlewRatesShortFoV_Filename,args.SlewRatesDiagonalFoV_Filename,args.SlewRatesLongFoV_Filename)

#Compute the yield for the input fields to test that everything is working.
if args.testyield or args.testplot or allfixed:
    totalYield, totalAreaPix, totalArea = handler.computeYield()
    print("Testing the yield computation.")
    print("Total yield:",totalYield)
    print("Total area (map pixels):",totalAreaPix)
    print("Total area (deg^2):",totalArea)

    if testplot:
        plt.figure()
        handler.plotFields()
        ym.plotMap()
        plt.title('Yield = %g',totalYield)
        plt.gca().set_aspect('equal')
        #if testplot==True:
        plt.show()
        #else:
        #    try:
        #        plt.savefig(testplot)
        #    except RuntimeError:
        #        raise('Error saving test figure')


#Store the originals
Cadence0 = args.Cadence0 + 0
texp0 = args.texp0 + 0
pathOverhead=0

if allfixed or allfree:
    #Compute the texp and cadence of the field layout once, no need
    #to recompute
    #except NotImplementedError:
    #    raise('Special case handling for allfixed or allfree not ')
    print("allfixed or allfree")
    handler.fromCentersChips(fields,romanFoV,ym)
    pathOverhead,bestPath = slewopt.optimizePath(fields)

    
if allfixed:
    cadence = (nread*args.readTime*Nfields+pathOverhead)/60.0
    handler.scaleMap(cadence,args.Cadence0,alphaC,
                     nread*args.readTime,args.texp0,alphaT)
    totalYield, totalAreaPix, totalArea = handler.computeYield()
    print("Total yield:",totalYield)
    print("Total area (map pixels):",totalAreaPix)
    print("Total area (deg^2):",totalArea)
    sys.exit()


#Build the grid to accept the results
if lrange[0]>lrange[1] and lstep>0:
    lstep = -lstep
if brange[0]>brange[1] and bstep>0:
    bstep = -bstep
lsteps = np.arange(lrange[0],lrange[1],lstep)
bsteps = np.arange(brange[0],brange[1],bstep)
lgrid,bgrid = np.meshgrid(lsteps,bsteps)
yieldgrid = np.empty(shape=lgrid.shape)
cadencegrid = np.empty(shape=lgrid.shape)
nreadgrid = np.empty(shape=lgrid.shape)

lcenter = np.average(fields['l'])
bcenter = np.average(fields['b'])
Nfields = fields.shape[0]

allBestFields = 0
allBestYield = -1e50

testfile = open(args.output_root + '_results.txt','w',buffering=1)


for index,l in np.ndenumerate(lgrid):

    b= bgrid[index]
    fieldsNew = fields.copy(deep=True)
    #print("fieldsNew before")
    #print(fieldsNew)
    fieldsNew.loc[freeMask,'l'] += (l-lcenter)
    fieldsNew.loc[freeMask,'b'] += (b-bcenter)
    #print("after")
    #print(fieldsNew)
    #print("")
    handler.fromCentersChips(fieldsNew,romanFoV,ym)
    if not allfree:
        #Recompute the path overhead if the distance between fields has
        #changed
        pathOverhead,bestPath = slewopt.optimizePath(fieldsNew)
    bestYield=-1e50
    bestCadence = np.nan
    bestNread = np.nan
    for nread in range(args.NRead_Bounds[0],args.NRead_Bounds[1]+1):
        cadence = (nread*args.readTime*Nfields+pathOverhead)/60.0
        #print(cadence)
        if args.Cadence_Bounds[0]<=cadence<args.Cadence_Bounds[1]:
            #Modify yieldMap for nread,cadence
            handler.scaleMap(cadence,args.Cadence0,alphaC,
                             nread*args.readTime,args.texp0,alphaT)
            #Compute the yield
            totalYield, totalAreaPix, totalArea = handler.computeYield()
            if totalYield > bestYield:
                bestCadence = cadence
                bestNread = nread
                bestYield = totalYield

            if totalYield > allBestYield:
                allBestYield = totalYield
                allBestFields = fieldsNew.copy(deep=True)
            #testfile.write("%g %g %d %g %g %g %g\n" % (l,b,nread,cadence,totalYield,totalAreaPix,totalArea))

    cadencegrid[index] = bestCadence
    nreadgrid[index] = bestNread
    yieldgrid[index] = bestYield
    print(l,b,bestNread*args.readTime,bestCadence,bestYield)
    testfile.write("%g %g %d %g %g\n" % (l,b,bestNread,bestCadence,bestYield))


handler.fromCentersChips(allBestFields,romanFoV,ym)
print("Best yield: ",allBestYield)
print("Best fields:")
print(allBestFields)
    
testfile.close()
with open(args.output_root + "_results.pkl",'wb') as pklhandle:
    pickle.dump([lgrid,bgrid,nreadgrid,cadencegrid,yieldgrid,handler],pklhandle)


        
    


