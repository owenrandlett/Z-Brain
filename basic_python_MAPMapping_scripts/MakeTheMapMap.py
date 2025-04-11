#%%

# This script will create the MAP-Map whole-brain activity map from
# registered pERK/tERK data by performing differential intensity analysis at
# each voxel. This is a python version of the original MakeTheMapMap.m script
# from Randlett et al 2015, Nature Methods. Results from this script are
# approximately the same as the original matlab script. 
#
# Morphed pERK data from groups to be compared needs to be in folders within a
# parent directory defined in 'root_dir', we are propted for this folder when
# the script starts (if we dont define it manually). Each parwise comparison
# between all the groups is made. The input images need to be single channel
# 8-bit tiff stacks. Any registered stacks registered to any reference brain
# should work, but to keep consistent with the formatting of the Z-Brain and
# subsequent analysis, images should be registered into the reference brain
# ('Ref20131120pt14pl2.nrrd'), downsampled to 300x679x80 (x,y,z), and
# pre-smoothed with a 2D gaussian filter (sigma = 2). This can be done in ImageJ
# with 'PrepareStacksForMAPMapping.ijm' or pytyon with
# 'PrepareStacksForMAPMapping.py' , which converts the '.nrrd' files output from
# CMTK or ANTS into Tiffs, downsample and smooths appropriately.
#
# We look for which files are pERK stains, and which are total ERK stains based
# on the naming, where channel 1 has the string '_01_' within it, while channel
# 2 is '_02_', etc. This needs to be set appropriately before running the
# script.

# Pairwise comparisons are made between all of the input folders. the two groups
# analyzed (GroupA and GroupB) are combined and then ranked. A Mann-Whitney U
# statistic Z score is calculated at each pixel. To then set an FDR threshold,
# pseudogroups are assembled by randomization such that they have an equal
# number of fish from GroupA and GroupB, randomly assigned, and a Z score is
# calculated as before. This mixing procedure is repeated 500 times, creating a
# distribution of control Zs at each voxel. From this distribution we then set a
# FDR threshold 0.00005 (0.005 probability of finding a value greater than this
# threshold). This is done assuming a normal distribution. The GroupA vs GroupB
# Z score is then thresholded at this value of Z (all pixels with lesser value
# are set to 0. The image is normalized such that 10 = 255 in an 8bit tiff, and
# the image stacks are written. I also set a requirement that there is a minimal
# pERK staining intensity at the pixel to get rid of most non-brain pixels
# (average across all brains must be pERK_thresh/255), and there also needs to be a
# mininmum of 5 warped fish in each group at the relevant pixel, set based on
# the fact that the only '0' pixels in the stack are those missing after warping
# - ie. not real image values.
#
#  The main file output from this analysis is the '*SignificantDeltaMedian.tif'
# file, which we refer to as the MAP-Map .  This stack depicts the delta Median
# value, at each pixel foud to be above  the FDR threshold. Significant pixels
# are assigned a non-0 value. Green signal indicates higher signals in the first
# Group (for example, GroupX in a 'GroupX_over_GroupY' comparison), magneta
# signals indicate higher signals  in the second group. These signals are mapped
# linearly between 0 and 65535 in a 16-bit stack, where 65535 represents a delta
# pERK level of 0.5.  We also write the thresholded Z-score comparison
# ('*ZScores.tif'), which represents is mapped 0-255 = ZScore 0-10 in an 8bit
# stack, As a control we write the final mixed group analysis from which the FDR
# threshold is calculated in order to give a sense of what false signals look
# like (there should be almost none). Finally we also write stacks depicting the
# median pERK level at each voxel, and the standard deviation at each voxel, for
# each group.
#
#--------------- Variables that can be tweaked
#
# terk_label and perk_label -- '*_0N_*' where N = post-warp channel number
#
# n_permutes -- the number of times we make pseudogroup comparisons for the 
#
#
# FDR calculation
#
# fdr_thresh -- The FRD thresholding of the Z-score. A value of 0.00005 gives
# very clean control comparisons.
#
# using_erk -- logical, set to True to run the total-ERK based normalization to
# calculate pERK levels (pERK/tERK). If you are not using this normalization
# (for example if you dont have the tERK stain, or want to quantify differences
# in another channel that should not be normalized, set to 0).
#
# write_mixed_control -- logical for whether to write the mixed control
# comparisons for a sanity check on the FDR threshold level Written by Owen
# Randlett (owen.randlett@gmail.com)

import numpy as np
import nrrd
import os
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import zoom, gaussian_filter
from scipy.stats import rankdata, norm
from skimage.io import imsave, imread
from natsort import natsorted
from PIL import Image
from tqdm.notebook import tqdm
from tkinter import Tk
from tkinter import filedialog
from multiprocessing import Pool, cpu_count
import time


# variables

terk_label = '/*_01_*'
perk_label = '/*_03_*'
n_permutes = 500
fdr_thresh = 0.00005 
pERK_thresh = 10 # minimum staining intensity to be considered
n_min_group = 5 # minumum number of fish in each group for a pixel to be condisered
n_cores = 25
using_erk = True
write_mixed_control = True

# separate folders for each group within one partent folder. Path of parent folder entered below as root directory, if left blank, will prompt for directory
root_dir = '/media/BigBoy/ciqle/LeicaStellaris/20250225_HuCGCAMP_pERK_tERK_Ethanol_5DF/SmoothedTiffs'  


if root_dir == '':
    root = Tk()
    root.withdraw()
    root_dir = filedialog.askdirectory(title='pick root directory containing folders with smoothed tiffs')

dir_out = os.path.join(root_dir, 'output_FDR='+str(fdr_thresh)+ '_UsingERK=' + str(using_erk)+ '/')
if not os.path.exists(dir_out):
    os.mkdir(dir_out)
#%
# functions

def get_mwu(ranks, ind_a, ind_b):

        # sum the ranks for each group
        ranksum_a = np.sum(ranks[ind_a, :, :], axis=0)
        ranksum_b = np.sum(ranks[ind_b, :, :], axis=0)

        #% determine the number of fish contributing to each pixel. This is not always the same because if for example part of the forebrain was not imaged in one fish, after warping there would be a 0-value pixels at this missing region. These are NaN'd out above, but we need to keep track of the Ns for the ranksum and z-score statistics
        n_in_a = len(ind_a) - np.sum(ranks[ind_a, :, :] == 0, axis=0)
        n_in_b = len(ind_b) - np.sum(ranks[ind_b, :, :] == 0, axis=0)

        # calculate the Mann-Whitney U-statistic (https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test)

        u_a = ranksum_a - (n_in_a * (n_in_a + 1))/2
        u_b = ranksum_b - (n_in_b * (n_in_b + 1))/2

        # calculate the Z-score. The distribution of Us can be assumed to be normal, and from inspection of my data this assumption appears valid
        # the value of U is the mean between the two groups
        mu = (u_a + u_b) / 2
        sig_u = np.sqrt( n_in_a*n_in_b * (n_in_a + n_in_b + 1) / 12 )
        zscores_a = (u_a - mu)/sig_u
        
        return zscores_a, n_in_a, n_in_b

def make_gm_img(save_name, im, sat_val):
    im_toconv = np.copy(im)
    im_toconv[~np.isfinite(im_toconv)] = 0
    im_toconv = im_toconv/sat_val
    im_toconv[im_toconv > 1] = 1
    im_toconv[im_toconv < -1] = -1
    im_pos = np.copy(im_toconv)
    im_pos[im_pos < 0] = 0
    im_neg = -np.copy(im_toconv)
    im_neg[im_neg < 0] = 0
    im_col = np.stack((im_neg, im_pos, im_neg), axis=3) * 65535
    im_col = im_col.astype('uint16')   
    imsave(save_name, im_col)


def compare_groups(compare):
    gr_a,dir_a, gr_b, dir_b = compare
    print(gr_a + ' vs ' + gr_b)
    # find the relevant files
    terk_names_a = natsorted(glob.glob(dir_a+terk_label))
    perk_names_a = natsorted(glob.glob(dir_a+perk_label))

    if not len(terk_names_a) == len(perk_names_a):
        raise ValueError('The number of stacks in group A arent equal')
    n_a = len(terk_names_a)

    terk_names_b = natsorted(glob.glob(dir_b+terk_label))
    perk_names_b = natsorted(glob.glob(dir_b+perk_label))

    if not len(terk_names_b) == len(perk_names_b):
        raise ValueError('The number of stacks in group B arent equal')

    n_b = len(terk_names_b)

    if n_a == 0: 
        raise ValueError(f'No groups found in directory {dir_a}')
    if n_b == 0: 
        raise ValueError(f'No groups found in directory {dir_b}')

    # preallocate images*
    im = imread(terk_names_a[0])
    im_dims = im.shape
    zs, height, width = im_dims

    # we will keep python convention of images with Zs, height, width as dimension order

    im_a_med = np.zeros(im_dims, dtype=float)
    im_a_std = np.zeros(im_dims, dtype=float)
    im_b_med = np.zeros(im_dims, dtype=float)
    im_b_std = np.zeros(im_dims, dtype=float)

    im_mwuz = np.zeros(im_dims, dtype=float)
    im_mwuz_singleit = np.zeros(im_dims, dtype=float)

    # set up the groups for the mixed group analysis

    rand_mix_a = np.zeros((int(np.floor(n_a / 2) * 2), n_permutes), dtype=int)
    rand_mix_b = np.zeros((int(np.floor(n_b / 2) * 2), n_permutes), dtype=int)

    for p in range(n_permutes):
        rand_mix_a[:, p] = np.random.choice(n_a, size=rand_mix_a.shape[0], replace=True)
        rand_mix_b[:, p] = n_a + np.random.choice(n_b, size=rand_mix_b.shape[0], replace=True)
    #%


    
    for z in tqdm(range(zs), desc=gr_a + ' vs ' + gr_b):
        # for each z-plane, we load in the relevant data and calculate the 

        # preallocate arrays to hold pERK data, and pERK ratio data
        im_a = np.zeros((n_a, height, width))
        im_a_perk = np.zeros((n_a, height, width))
        im_b = np.zeros((n_b, height, width))
        im_b_perk = np.zeros((n_b, height, width))

        # load the relevant slices
        for n in range(n_a):
            im_a_perk[n,:,:] = imread(perk_names_a[n], img_num=z)
            if using_erk:
                im_a[n,:,:] = im_a_perk[n,:,:] / imread(terk_names_a[n], img_num=z)
            else:
                im_a[n,:,:] = im_a_perk[n,:,:]

        for n in range(n_b):
            im_b_perk[n,:,:] = imread(perk_names_b[n], img_num=z)
            if using_erk:
                im_b[n,:,:] = im_b_perk[n,:,:] / imread(terk_names_b[n], img_num=z)
            else:
                im_b[n,:,:] = im_b_perk[n,:,:]

        im_a_med[z, :, :] = np.nanmedian(im_a, axis=0)
        im_b_med[z, :, :] = np.nanmedian(im_b, axis=0)

        im_a_std[z, :, :] = np.nanstd(im_a, axis=0)
        im_b_std[z, :, :] = np.nanstd(im_b, axis=0)

        im_comb = np.vstack((im_a, im_b))

        # there are NaNs and Inf values in the im_comb array due to 0 values in the images. We will assume that there instaces of 0 are areas of the fish that were not imaged and so we will ignore them from the analysis by setting to 0.
        im_comb[~np.isfinite(im_comb)] = np.nan
        im_comb[im_comb == 0] = np.nan

        ranks = rankdata(im_comb, axis=0, method='average', nan_policy = 'omit')
        
        ranks[~np.isfinite(im_comb)] = 0

        zscores_mixed_a = np.zeros((n_permutes, height, width))
        zscores_mixed_b = np.zeros((n_permutes, height, width))

        

        for p in range(n_permutes):
            
            # get the vector for the indieices of the fish in this mixed group for this permutation
            rand_a_perm = rand_mix_a[:,p] 
            rand_b_perm = rand_mix_b[:,p]

            fish_mixed_a = np.hstack((
                rand_a_perm[:int(len(rand_a_perm)/2)], 
                rand_b_perm[:int(len(rand_b_perm)/2)]))

            fish_mixed_b = np.hstack((
                rand_a_perm[int(len(rand_a_perm)/2):], 
                rand_b_perm[int(len(rand_b_perm)/2):]))
            
            zscores_mixed_a[p,:,:], n_in_a_mix, n_in_b_mix =  get_mwu(ranks, fish_mixed_a, fish_mixed_b)


        # now we set the FDR threshold. The distribution of Zs appears normal at each pixel, so I will use the inverse normal cdf

        zscores_thresh = norm.ppf(1-fdr_thresh, loc=0, scale=np.std(zscores_mixed_a, axis=0))

        # now run the actual comparison between the two groups
        fish_a = np.arange(n_a)
        fish_b = np.arange(n_a, n_a+n_b)
        zscores_a, n_in_a, n_in_b = get_mwu(ranks, fish_a, fish_b)

        # throw away any pixel where the mean pERK value is less than the pERK threshold, default = 10 (as the brain is highly stained, these are mostly non-brain pixels)

        mean_perk = np.nanmean(np.vstack((im_a_perk, im_b_perk)), axis=0)
        zscores_a[mean_perk < pERK_thresh] = np.nan

        # throw away any pixel where we dont have at least n_min_group fish from both groups represented
        zscores_a[n_in_a < n_min_group] = np.nan
        zscores_a[n_in_b < n_min_group] = np.nan

        # threshold the z-scores based on the FDR threshold we calculated, note that this is done separately on each pixel

        zscores_a[abs(zscores_a)<zscores_thresh] = np.nan

        # assign these values in the image stack we will write later
        im_mwuz[z, :, :] = zscores_a

        # do the same thing for the last mixed-group comparison from the bootstrapping as an example of the false discoveries

        zscores_singleit = zscores_mixed_a[-1,:,:]
        zscores_singleit[mean_perk < pERK_thresh] = np.nan
        zscores_singleit[n_in_a_mix < n_min_group] = np.nan
        zscores_singleit[n_in_b_mix < n_min_group] = np.nan
        zscores_singleit[abs(zscores_singleit)<zscores_thresh] = np.nan

        im_mwuz_singleit[z,:,:] = zscores_singleit
       

    os.chdir(dir_out)
    make_gm_img(f'{gr_a}_over_{gr_b}_ZScores.tif', im_mwuz, 10)
    make_gm_img(f'{gr_b}_over_{gr_a}_ZScores.tif', -im_mwuz, 10)

    # compulte delta median, thresholded with z-scores above FDR threshold
    delta_med = (im_a_med - im_b_med)/(im_a_med + im_b_med)
    delta_med[~np.isfinite(im_mwuz)] = 0

    make_gm_img(f'{gr_a}_over_{gr_b}_SignificantDeltaMedians.tif', delta_med, 2)
    make_gm_img(f'{gr_b}_over_{gr_a}_SignificantDeltaMedians.tif', -delta_med, 2)

    if write_mixed_control:
        make_gm_img(f'MixedControl_{gr_a}_over_{gr_b}_ZScores.tif', im_mwuz_singleit, 10)
        delta_med_mixed = (im_a_med - im_b_med)/(im_a_med + im_b_med)
        delta_med_mixed[~np.isfinite(im_mwuz_singleit)] = 0
        make_gm_img(f'MixedControl_{gr_a}_over_{gr_b}_SignificantDeltaMedians.tif', delta_med_mixed, 2)

    # write the median stack
    if using_erk: # write as 16 bit images normalized to 10
        imsave(f'{gr_a}_MedianImage.tif', (65535*im_a_med/10).astype('uint16'))
        imsave(f'{gr_b}_MedianImage.tif', (65535*im_b_med/10).astype('uint16'))
        imsave(f'{gr_a}_StdImage.tif', (65535*im_a_std/10).astype('uint16'))
        imsave(f'{gr_b}_StdImage.tif', (65535*im_b_std/10).astype('uint16'))
    else: # if not using ERK these will be 8 bit images
        imsave(f'{gr_a}_MedianImage.tif', im_a_med.astype('uint8'))
        imsave(f'{gr_b}_MedianImage.tif', im_b_med.astype('uint8'))
        imsave(f'{gr_a}_StdImage.tif', im_a_std.astype('uint8'))
        imsave(f'{gr_b}_StdImage.tif', im_b_std.astype('uint8'))


# sort out directories and perform comparisons


sub_dirs = os.listdir(root_dir)
dirs = []
names = []
for i, sub_dir in enumerate(sub_dirs):
    if sub_dir.find('output') == -1 and sub_dir.find('Exclude'): # ignore output directories and "Exclude" directory, if they exist
        dirs.append(os.path.join(root_dir,sub_dir))
        names.append(sub_dir)
        print('Group was detected as ' + sub_dir)

n_groups = len(names)
n_compare = int(n_groups*(n_groups-1)/2)#
comparisons = []
for m in np.arange(n_groups):
    for n in np.arange(m+1, n_groups):
        gr_a = names[m]
        dir_a = dirs[m]
        gr_b = names[n]
        dir_b =  dirs[n]
        comparisons.append([gr_a, dir_a, gr_b, dir_b])

print('... starting comparisons ...')

#%
n_cores = np.min((n_compare, n_cores))
with Pool(min(n_cores, cpu_count())) as p:
    p.map(compare_groups, comparisons)
# for compare in comparisons:
#     compare_groups(compare)

#%%

