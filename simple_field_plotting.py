import matplotlib.pyplot as plt
import pandas as pd

sca = pd.read_csv('outline_sca_layout.txt',sep='\s+',header=None)
fields = pd.read_csv('field_layouts/layout_163000.centers',sep='\s+')

plt.figure()

for i,f in fields.iterrows():
    ls='k-'
    if f['fixed']==1:
        ls='r-'

    plt.plot(sca.iloc[:,1]+f['l'],sca.iloc[:,2]+f['b'],ls)
    

plt.xlabel('l (deg)')
plt.ylabel('b (deg)')
plt.gca().invert_xaxis()
plt.gca().set_aspect('equal')
plt.show()
