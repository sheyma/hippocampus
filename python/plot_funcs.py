import os
import h5py
import numpy as np
from scipy.io import loadmat

# get surface infos

ddir = '../data/tout_group/'

# get surface infos
LSUBfile = os.path.join(ddir, 'surf_lsub.mat')
LSUBsurf = loadmat(LSUBfile)['ave_lsub'][0,0]
xLSUB = LSUBsurf['coord'][0,:]           # (1024,)
yLSUB = LSUBsurf['coord'][1,:]           # (1024,)
zLSUB = -LSUBsurf['coord'][2,:]          # (1024,)
triLSUB = LSUBsurf['tri']                 # (2044, 3) 

LCAfile = os.path.join(ddir, 'surf_lca.mat')
LCAsurf = loadmat(LCAfile)['ave_lca'][0,0];
xLCA = LCAsurf['coord'][0,:]           # (2048,)
yLCA = LCAsurf['coord'][1,:]           # (2048,)
zLCA = -LCAsurf['coord'][2,:]          # (2048,)
triLCA = LCAsurf['tri']                 # (4092, 3)

LDGfile = os.path.join(ddir, 'surf_ldg.mat')
LDGsurf = loadmat(LDGfile)['ave_ldg'][0,0];
xLDG = LDGsurf['coord'][0,:]           # (1024,)
yLDG = LDGsurf['coord'][1,:]           # (1024,)
zLDG = -LDGsurf['coord'][2,:]          # (1024,)
triLDG = LDGsurf['tri']                 # (2044, 3)

# get surface infos
RSUBfile = os.path.join(ddir, 'surf_rsub.mat')
RSUBsurf = loadmat(RSUBfile)['ave_lsub'][0,0]
xRSUB = RSUBsurf['coord'][0,:]           # (1024,)
yRSUB = RSUBsurf['coord'][1,:]           # (1024,)
zRSUB = -RSUBsurf['coord'][2,:]          # (1024,)
triRSUB = RSUBsurf['tri']                 # (2044, 3) 

RCAfile = os.path.join(ddir, 'surf_rca.mat')
RCAsurf = loadmat(RCAfile)['ave_lca'][0,0];
xRCA = RCAsurf['coord'][0,:]           # (2048,)
yRCA = RCAsurf['coord'][1,:]           # (2048,)
zRCA = -RCAsurf['coord'][2,:]          # (2048,)
triRCA = RCAsurf['tri']                 # (4092, 3)

RDGfile = os.path.join(ddir, 'surf_rdg.mat')
RDGsurf = loadmat(RDGfile)['ave_ldg'][0,0];
xRDG = RDGsurf['coord'][0,:]           # (1024,)
yRDG = RDGsurf['coord'][1,:]           # (1024,)
zRDG = -RDGsurf['coord'][2,:]          # (1024,)
triRDG = RDGsurf['tri']                 # (2044, 3)



def plot_surf(xS, yS, zS, triS, data, cmap, cmin, cmax):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import seaborn as sns
    sns.set_context("poster", font_scale=1.4)
    sns.set_style('white')
    
    x = xS          
    y = yS
    z = zS       

    triangles = triS - 1 # matlab-python index matching 

    c = data

    colors = np.mean( [c[triangles[:,0]], c[triangles[:,1]], c[triangles[:,2]]], axis = 0);

    # Displays the 4D graphic.
    fig = plt.figure(figsize=(12,12));
    ax  = fig.add_subplot(1, 1, 1, projection='3d')
    ax.view_init(azim=240)

    surf = ax.plot_trisurf(x,y,z,
                           triangles = triangles, 
                           cmap = cmap,
                           #cmap = 'Greys_r',
                           #cmap = parula_map,
                           edgecolor = None,
                           linewidth = 0, 
                           antialiased = False)

    surf.set_array(colors)
    surf.autoscale()
    surf.set_clim(cmin, cmax)
    
    cbar = fig.colorbar(surf, shrink=0.5, aspect=10, ticks=[0.5, 0.8, 1.1])
    cbar.ax.get_yaxis().labelpad = 15 
    cbar.ax.tick_params(labelsize = 60)
    ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_axis_off()
    
    return fig




def plot_surf_upper2(xS, yS, zS, triS, data, cmap, cmin, cmax):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import seaborn as sns

    g = sns.set_style('white')
    
    x = xS          
    y = yS
    z = zS       

    triangles = triS - 1 # matlab-python index matching 

    c = data

    colors = np.mean( [c[triangles[:,0]], c[triangles[:,1]], c[triangles[:,2]]], axis = 0);

    # Displays the 4D graphic.
    fig = plt.figure(figsize=(12,12));
    ax  = fig.add_subplot(1, 1, 1, projection='3d')
    ax.view_init(azim=240)

    surf = ax.plot_trisurf(x,y,z,
                           triangles = triangles, 
                           cmap = cmap,
                           #cmap = 'Greys_r',
                           #cmap = parula_map,
                           edgecolor = None,
                           linewidth = 0, 
                           antialiased = False)

    surf.set_array(colors)
    surf.autoscale()
    surf.set_clim(cmin, cmax)
    
    cbar = fig.colorbar(surf, shrink=0.5, aspect=10, ticks=[cmin,  cmax])
    cbar.ax.set_yticklabels(["{:1.2f}".format(i) for i in [cmin, cmax]])

    cbar.ax.get_yaxis().labelpad = 15 
    cbar.ax.tick_params(labelsize = 60)    
    
    ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_axis_off()
    
    l, b, w, h = ax.get_position().bounds
    ll, bb, ww, hh = cbar.ax.get_position().bounds
    ax.set_position([l-0.05, b, w, h])
    cbar.ax.set_position([ll-0.1, b + 0.1*h, ww, h*0.8])

    
    return fig


def plot_surf_upper3(xS, yS, zS, triS, data, cmap, cmin, cmax):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import seaborn as sns

    g = sns.set_style('white')
    
    x = xS          
    y = yS
    z = zS       

    triangles = triS - 1 # matlab-python index matching 

    c = data

    colors = np.mean( [c[triangles[:,0]], c[triangles[:,1]], c[triangles[:,2]]], axis = 0);

    # Displays the 4D graphic.
    fig = plt.figure(figsize=(12,12));
    ax  = fig.add_subplot(1, 1, 1, projection='3d')
    ax.view_init(azim=240)

    surf = ax.plot_trisurf(x,y,z,
                           triangles = triangles, 
                           cmap = cmap,
                           #cmap = 'Greys_r',
                           #cmap = parula_map,
                           edgecolor = None,
                           linewidth = 0, 
                           antialiased = False)

    surf.set_array(colors)
    surf.autoscale()
    surf.set_clim(cmin, cmax)
    
    cbar = fig.colorbar(surf, shrink=0.5, aspect=10, ticks=[cmin,  cmax])
    cbar.ax.set_yticklabels(["{:1.3f}".format(i) for i in [cmin, cmax]])

    cbar.ax.get_yaxis().labelpad = 15 
    cbar.ax.tick_params(labelsize = 60)    
    
    ax.grid(False)

    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_axis_off()
    
    l, b, w, h = ax.get_position().bounds
    ll, bb, ww, hh = cbar.ax.get_position().bounds
    ax.set_position([l-0.05, b, w, h])
    cbar.ax.set_position([ll-0.1, b + 0.1*h, ww, h*0.8])

    
    return fig




def parula_cmap():

    from matplotlib.colors import LinearSegmentedColormap
    
    cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
     [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
     [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
      0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
     [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
      0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
     [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
      0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
     [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
      0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
     [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
      0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
     [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
      0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
      0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
     [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
      0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
     [0.0589714286, 0.6837571429, 0.7253857143], 
     [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
     [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
      0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
     [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
      0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
     [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
      0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
     [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
      0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
     [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
     [0.7184095238, 0.7411333333, 0.3904761905], 
     [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
      0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
     [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
     [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
      0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
     [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
      0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
     [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
     [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
     [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
      0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
     [0.9763, 0.9831, 0.0538]]
    
    parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
    
    return parula_map
    

    
import matplotlib.colors as colors
import matplotlib.pyplot as plt

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = plt.get_cmap('nipy_spectral')
new_cmap = truncate_colormap(cmap, 0.2, 0.95)

#colors1 = plt.cm.YlGnBu(np.linspace(0, 1, 128))
first = 70
second = 58
third  = 64
fourth = 64

colors1 = plt.cm.jet(np.linspace(0.67, 0.7, third))
colors2 = plt.cm.hot(np.linspace(0, 0.8, fourth))

cols = np.vstack((colors2))
hotcolors = colors.LinearSegmentedColormap.from_list('my_colormap', cols)

#num = 256
#gradient = range(num)
#for x in range(5):
#    gradient = np.vstack((gradient, gradient))

#fig, ax = plt.subplots(nrows=1)
#ax.imshow(gradient, cmap=pnasEDITED, interpolation='nearest')
#ax.set_axis_off()
#fig.tight_layout()

#plt.show()


