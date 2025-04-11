#%% this script makes projections of the MAPMAP output files in a 2d image of the Z and X projections

from skimage.io import imsave, imread
from scipy.ndimage import zoom
import os, glob, natsort
import numpy as np
import matplotlib.pyplot as plt


dir_mapmap = r'/media/BigBoy/ciqle/LeicaStellaris/20250225_HuCGCAMP_pERK_tERK_Ethanol_5DF/SmoothedTiffs/output_FDR=0.0005_UsingERK=True'
file_ZBrain2_0_outline_proj = r'/home/zeneb/github/fish_registration/basic_MAPMapping_Scripts/zbrain2.0/ZBrain2_0_outline_proj.tif'

dir_mapmap = os.path.realpath(dir_mapmap)

calc_max_val = True
max_val = 1250

#dir_proj = os.path.realpath(r'/media/BigBoy/Owen/GranatoLabData/20180113_pERK_Hab_Melatonin_Picro_Hexestrol/SmoothedTiffs/output_FDR=5e-05')
# make the projections: 


dir_out = os.path.join(dir_mapmap, 'projections/')
if not os.path.exists(dir_out):
    os.mkdir(dir_out)


orig_rez = np.array([2, 0.798, 0.798])
orig_size =  np.array([138, 1406,  621])
IM_outlines = imread(file_ZBrain2_0_outline_proj)

os.chdir(dir_mapmap)
comp_files = glob.glob('*_SignificantDeltaMedians.tif')
comp_files_noMixed = []
for file in comp_files:
    if file.find('MixedControl') == -1:
        comp_files_noMixed.append(file)
        print(file)

for file in comp_files_noMixed:
    IM  = imread(file).astype('double')
    new_size = IM.shape[:3]

    scale_factor = orig_size/new_size
    new_rez = orig_rez * scale_factor

    IM_z = np.mean(IM, axis=0)
    IM_x = np.moveaxis(np.mean(IM, axis=2), 1, 0)
    IM_x_zoom = zoom(IM_x, [1,new_rez[0]/new_rez[2], 1], order =1)


    IM_proj = np.hstack((IM_z,IM_x_zoom))
    if calc_max_val:
        max_val = np.max(IM_proj)

    IM_proj = IM_proj/max_val
    IM_outlines_zoom = zoom(IM_outlines, np.array(IM_proj.shape[:2])/ np.array(IM_outlines.shape), order=0)
    imsave(file.replace('.tif', '_projection_maxval-' + str(int(max_val)) +'_.tif'), IM_proj)
    for ch in range(3):
        sl = IM_proj[:,:,ch]
        sl[IM_outlines_zoom > 0] = 1
        IM_proj[:,:,ch] = sl
    imsave(os.path.join(dir_out, file.replace('.tif', '_projection_maxval-' + str(int(max_val)) +'_wOutlines.tif')), IM_proj)
    plt.figure(figsize=(10,30))
    plt.imshow(IM_proj*5)
    plt.title(file)
    plt.show()

print('done')

#%% make a collage of certain subsets
comp_files_noMixed = natsort.natsorted(comp_files_noMixed)
for k, file in enumerate(comp_files_noMixed):
    print('...')
    print(k)
    print(file)

#%%

inds = [11, 7, 9]

comp_files_subset=[]
for ind in inds:
    comp_files_subset.append(comp_files_noMixed[ind])
comp_title = ''

print('..............................')
print('selected files = ')

print(comp_files_subset)

for k, file in enumerate(comp_files_subset):
    comp_title = comp_title + file.replace('SignificantDeltaMedians.tif', ':--:')
    IM  = imread(file).astype('double')
    new_size = IM.shape[:3]

    scale_factor = orig_size/new_size
    new_rez = orig_rez * scale_factor

    IM_z = np.mean(IM, axis=0)
    IM_x = np.moveaxis(np.mean(IM, axis=2), 1, 0)
    IM_x_zoom = zoom(IM_x, [1,new_rez[0]/new_rez[2], 1], order =1)


    IM_proj = np.hstack((IM_z,IM_x_zoom))
    max_val = np.max(IM_proj)
    #IM_proj = IM_proj/max_val
    IM_outlines_zoom = zoom(IM_outlines, np.array(IM_proj.shape[:2])/ np.array(IM_outlines.shape), order=0)
    

    for ch in range(3):
        sl = IM_proj[:,:,ch]
        sl[IM_outlines_zoom > 0] = 65535
        IM_proj[:,:,ch] = sl
    

    if k == 0:
        IM_proj_comb = IM_proj
        max_val_total = max_val
    else:
        IM_proj_comb = np.hstack((IM_proj_comb, IM_proj))
        max_val_total = np.max((max_val_total, max_val))

IM_proj_comb = IM_proj_comb/max_val_total
IM_proj_comb[IM_proj_comb > 1] = 1

plt.figure(figsize=(30,30))
plt.imshow(IM_proj_comb*3)
plt.title(comp_title)
plt.show()


imsave(comp_title + '_projection_maxval-' + str(int(max_val_total)) +'_wOutlines.tif', IM_proj_comb)
# %%
