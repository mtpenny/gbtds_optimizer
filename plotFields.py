'''A minimal working example for plotting fields with chips'''


from fields import fov, fovHandler
import argparse
from yieldMap import yieldMap
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Optimizer for the Roman Galactic Bulge Time Domain Survey")
parser.add_argument('yieldmap_filename',
                    help='File name for the yield map')
parser.add_argument('fields_filename',
                    help='File containing the field specifications (columns: name l b fixed)')
parser.add_argument('--sca-filename',default='sca_layout.txt',
                    help='Filename of the sca vertices')
args = parser.parse_args()


ym = yieldMap(args.yieldmap_filename)


#Load the field vertices
romanFoV = fov(args.sca_filename,unit='deg')

try:
    fields = pd.read_csv(args.fields_filename,sep='\s+',header=None,names=['field','l','b','fixed'],
                         dtype={'field':str,'l':float,'b':float,'fixed':int})
except:
    try:
        fields = pd.read_csv(args.fields_filename,sep='\s+')
    except:
        raise RuntimeError('Error reading fields file (%s) - file does not exist or is in an incorrect format. It should containt the named columns ["field","l","b","fixed"] or 4 unnamed columns' % (args.fields_filename))

handler = fovHandler()

handler.fromCentersChips(fields,romanFoV,ym)

plt.figure()
#Can pass this an ax if desired in a subplot
handler.plotFields(plot_kwargs={'color':'r','linestyle':'-','lw':0.5})
plt.show()
