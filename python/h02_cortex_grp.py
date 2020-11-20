"""
calculates the group-level connectivity from cortex to hippocampus
usage: $ python h02_cortex_grp.py
"""


import os
import h5py
import numpy as np
from numpy import genfromtxt    
from brainspace.datasets import load_conte69
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels

ddir      = '../data/'                                          #  data dir
cordir    = '../data/tout_cortex/'
odir      = '../data/tout_group'

# final subject list after QC         
subjlist = os.path.join(ddir, 'subjectListS900_QC_gr.txt');     # 709 subjects
f = open(subjlist); mylist = f.read().split("\n"); f.close() 
mylist =  mylist[:-1] 
totnum = len(mylist)

print('We have now %i subjects... ' % totnum)

avg_LSUB = np.zeros((360,1))
avg_LCA  = np.zeros((360,1))
avg_LDG  = np.zeros((360,1))
avg_RSUB = np.zeros((360,1))
avg_RCA  = np.zeros((360,1))
avg_RDG  = np.zeros((360,1))

for subjID in mylist:   
    # get the cortex-sub connectivity for each subject & append
    subjsub = os.path.join(cordir, subjID + '_cortex_LSUB.h5')
    with h5py.File(subjsub, "r") as f:        
        subjdata = np.array(f[subjID])   
        avg_LSUB = avg_LSUB + subjdata

   # get the cortex-ca connectivity for each subject & append
    subjca = os.path.join(cordir, subjID + '_cortex_LCA.h5')
    with h5py.File(subjca, "r") as f:        
        subjdata = np.array(f[subjID])   
        avg_LCA = avg_LCA + subjdata

    # get the cortex-dg connectivity for each subject & append
    subjdg = os.path.join(cordir, subjID + '_cortex_LDG.h5')
    with h5py.File(subjdg, "r") as f:        
        subjdata = np.array(f[subjID])   
        avg_LDG = avg_LDG + subjdata

    # get the cortex-sub connectivity for each subject & append
    subjsubR = os.path.join(cordir, subjID + '_cortex_RSUB.h5')
    with h5py.File(subjsubR, "r") as f:        
        subjdata = np.array(f[subjID])   
        avg_RSUB = avg_RSUB + subjdata

   # get the cortex-ca connectivity for each subject & append
    subjcaR = os.path.join(cordir, subjID + '_cortex_RCA.h5')
    with h5py.File(subjcaR, "r") as f:        
        subjdata = np.array(f[subjID])   
        avg_RCA = avg_RCA + subjdata

    # get the cortex-dg connectivity for each subject & append
    subjdgR = os.path.join(cordir, subjID + '_cortex_RDG.h5')
    with h5py.File(subjdgR, "r") as f:        
        subjdata = np.array(f[subjID])   
        avg_RDG = avg_RDG + subjdata

avg_LSUB = avg_LSUB / float(totnum)
avg_LCA  = avg_LCA / float(totnum)
avg_LDG  = avg_LDG / float(totnum)
avg_RSUB = avg_RSUB / float(totnum)
avg_RCA  = avg_RCA / float(totnum)
avg_RDG  = avg_RDG / float(totnum)

print(avg_LSUB.min(), avg_LSUB.max())
print(avg_LCA.min(), avg_LCA.max())
print(avg_LDG.min(), avg_LDG.max())
print(avg_RSUB.min(), avg_RSUB.max())
print(avg_RCA.min(), avg_RCA.max())
print(avg_RDG.min(), avg_RDG.max())

h = h5py.File(os.path.join(odir, 'cortex709_LSUB.h5'), 'w')
h.create_dataset('LSUB', data = avg_LSUB)
h.close()

h = h5py.File(os.path.join(odir, 'cortex709_LCA.h5'), 'w')
h.create_dataset('LCA', data = avg_LCA)
h.close()

h = h5py.File(os.path.join(odir, 'cortex709_LDG.h5'), 'w')
h.create_dataset('LDG', data = avg_LDG)
h.close()

h = h5py.File(os.path.join(odir, 'cortex709_RSUB.h5'), 'w')
h.create_dataset('RSUB', data = avg_RSUB)
h.close()

h = h5py.File(os.path.join(odir, 'cortex709_RCA.h5'), 'w')
h.create_dataset('RCA', data = avg_RCA)
h.close()

h = h5py.File(os.path.join(odir, 'cortex709_RDG.h5'), 'w')
h.create_dataset('RDG', data = avg_RDG)
h.close()

