'''A minimal working example for plotting fields with chips'''


from fields import fov, fovHandler
import argparse
from yieldMap import yieldMap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Optimizer for the Roman Galactic Bulge Time Domain Survey")
parser.add_argument('yieldmap_filename',
                    help='File name for the yield map')
parser.add_argument('fields_filename',
                    help='File containing the field specifications (columns: name l b fixed)')
parser.add_argument('--sca-filename',default='sca_layout.txt',
                    help='Filename of the sca vertices')
parser.add_argument('--location',nargs=2,default=[None,None],type=float,
                    help='Placement of the field - provide l and b as arguments. If no arguments, it is assumed the field will be at the position in the layout file')
parser.add_argument('--out-file',default='test.covfac',
                    help='Filename to store the covfac file')
parser.add_argument('--plot',action='store_true',
                    help='Plot the covfac')

args = parser.parse_args()


ym = yieldMap(args.yieldmap_filename)


#Load the field vertices
romanFoV = fov(args.sca_filename,unit='deg')

try:
    fields = pd.read_csv(args.fields_filename,sep=r'\s+',header=None,names=['field','l','b','fixed'],
                         dtype={'field':str,'l':float,'b':float,'fixed':int})
except:
    try:
        fields = pd.read_csv(args.fields_filename,sep=r'\s+')
    except:
        raise RuntimeError('Error reading fields file (%s) - file does not exist or is in an incorrect format. It should containt the named columns ["field","l","b","fixed"] or 4 unnamed columns' % (args.fields_filename))

if args.location[0] is not None:

    allfixed = False
    allfree = True
    #fixedMask = fields['fixed'].str.contains('fixed')
    fixedMask = fields['fixed']>0
    freeMask = np.logical_not(fixedMask)
    nfixed = fixedMask.sum()
    if nfixed==fields.shape[0]:
        allfixed = True
    if nfixed>0:
        allfree = False


    if allfixed:
        lcenter = 0.5*(np.max(fields['l'])+np.min(fields['l']))
        bcenter = 0.5*(np.max(fields['b'])+np.min(fields['b']))
    else:
        lcenter = 0.5*(np.max(fields.loc[freeMask,'l'])+np.min(fields.loc[freeMask,'l']))
        bcenter = 0.5*(np.max(fields.loc[freeMask,'b'])+np.min(fields.loc[freeMask,'b']))
    
    fields['l'] += args.location[0]-lcenter
    fields['b'] += args.location[1]-bcenter

handler = fovHandler()

handler.fromCentersChips(fields,romanFoV,ym)

handler.yieldMap.covfac = np.zeros(handler.yieldMap.yieldmap.shape)
totalYieldtmp, totalAreaPixtmp, totalAreatmp = handler.computeYield(True)

covfac = pd.DataFrame({'l':handler.yieldMap.llist,'b':handler.yieldMap.blist,'covfac':handler.yieldMap.covfac.flatten()})
covfac.to_csv(args.out_file,index=False)
                    


if args.plot:
    plt.figure()
    #Can pass this an ax if desired in a subplot
    plt.scatter(covfac['l'],covfac['b'],c=covfac['covfac'],s=20,marker='s')
    plt.gca().set_aspect('equal')
    plt.show()
