import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

sca_layout = pd.read_csv('sca_layout.txt',sep=r'\s+',header=None,skip_blank_lines=True,names=['n','l','b'])


extra_fields = pd.DataFrame({"name":["gal-center","globular-clusters","etz"],
                             "field":["GC","GLOB","ETZ"],
                             "l":[0.0,0.0,6.4],
                             "b":[-0.125,-4.2,0.0],
                             "fixed":[1,1,1]})

sca_row = sca_layout[(sca_layout['n']>=10) & (sca_layout['n']<=12)]
lwidth = sca_row['l'].max()-sca_row['l'].min()
bheight = sca_layout['b'].max()-sca_layout['b'].min()

print(lwidth,bheight)


start_point = [1.0,-1.0]
letters = 'abcdefghijklmnopqrstuvwxyz'
nfields = [4, 5, 6, 7, 8, 9]


for nf in nfields:
    try:
        os.mkdir(f'field_layouts/{nf}fields/')
    except:
        pass

    ids = np.arange(nf).astype(int)
    for nupper in range(nf):
        for let,ushift in enumerate(np.arange(-1.0,(ids[-1]-nupper)+1.1,0.5)*lwidth):
            l = np.zeros(ids.shape) + start_point[0]
            b = np.zeros(ids.shape) + start_point[1]
            if nupper>0:
                l[ids<nupper] += (ids*lwidth+ushift)[ids<nupper]
                b[ids<nupper] += bheight
            l[ids>=nupper] += ((ids-nupper)*lwidth)[ids>=nupper]

            fixed=np.zeros(ids.shape)

            output = pd.DataFrame()
            output['field'] = ids
            output['l'] = l
            output['b'] = b
            output['fixed'] = fixed

            #print(output)

            fname = f'field_layouts/{nf}fields/layout_{nf}f_{nupper}{letters[let]}.centers'
            print(fname)
            output.to_csv(fname,sep=' ')

            #print(pd.({"field":"","l":1000.0,"b":1000.0,"fixed":1}))

            output.loc[nf] = pd.Series({"field":"","l":1000.0,"b":1000.0,"fixed":1})
            for ef,efrow in extra_fields.iterrows():
                output.loc[nf] = efrow[["field","l","b","fixed"]]
                fname=fname = f'field_layouts/{nf}fields/layout_{nf}f_{nupper}{letters[let]}_{efrow["name"]}.centers'
                #print(f'field_layouts/{nf}fields/layout_{nf}f_{nupper}{letters[let]}_{efrow["name"]}.centers')
                #print(output)
                print(fname)
                output.to_csv(fname,sep=' ')

        if nupper==0:
            continue











