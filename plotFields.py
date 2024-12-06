'''A minimal working example for plotting fields with chips'''

import numpy as np
from fields import fov, fovHandler
import argparse
from yieldMap import yieldMap
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

parser = argparse.ArgumentParser(description="Optimizer for the Roman Galactic Bulge Time Domain Survey")
parser.add_argument('yieldmap_filename',
                    help='File name for the yield map')
parser.add_argument('fields_filename',
                    help='File containing the field specifications (columns: name l b fixed)')
parser.add_argument('--sca-filename',default='sca_layout.txt',
                    help='Filename of the sca vertices')
parser.add_argument('--save',default=None,
                    help='Filename of the sca vertices')
parser.add_argument('--location',nargs=2,default=[None,None],type=float,
                    help='Placement of the field - provide l and b as arguments. If no arguments, it is assumed the field will be at the position in the layout file')
parser.add_argument('--lrange',nargs=2,default=[3.0,-3.0],type=float,
                    help='Range of l')
parser.add_argument('--brange',nargs=2,default=[-3.0,3.0],type=float,
                    help='Range of b')
parser.add_argument('--cblabel',default=r'Yield per survey per map tile',type=str,
                    help='Change the text of the colorbar label')

args = parser.parse_args()

location = [0.0,0.0]
if not args.location[0] is None and not args.location[1] is None:
    location = args.location

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

print("Shift:",location[0],location[1])
    
notfixed = (fields['fixed'] == 0)
#print(fields[['l','b']])
#print(notfixed)
shiftl = np.zeros(fields['l'].shape)
shiftb = np.zeros(fields['b'].shape)
shiftl[notfixed] = location[0]
shiftb[notfixed] = location[1]
fields['l'] += shiftl
fields['b'] += shiftb
print("Fields:")
print(fields[['l','b','fixed']])

handler = fovHandler()

handler.fromCentersChips(fields,romanFoV,ym)

plt.figure()
#plt.rcParams["font.family"] = 'serif'
#plt.rcParams["font.serif"] = 'Times'
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams["font.family"] = 'STIXGeneral'
plt.rcParams["font.size"] = 16
#plt.rcParams['lines.linewidth'] = 5
lw = 3
plt.rcParams['axes.linewidth'] = lw
#Can pass this an ax if desired in a subplot
ymap = handler.yieldMap.plotMap(pcolormesh_kwargs={'cmap':cm.gray_r})
handler.plotFields(plot_kwargs={'color':'r','linestyle':'-','lw':2})
plt.xlim(args.lrange)
plt.ylim(args.brange)
plt.gca().set_aspect('equal')
plt.xlabel('l (deg)')
plt.ylabel('b (deg)')
plt.colorbar(label=args.cblabel)
ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax.tick_params(axis='x',direction='in',length=10,width=lw/2.0,bottom=True,top=True)
ax.tick_params(axis='y',direction='in',length=10,width=lw/2.0,left=True,right=True)

plt.tight_layout()


if args.save is None:
    plt.show()
else:
    print(args.save)
    plt.savefig(args.save)
