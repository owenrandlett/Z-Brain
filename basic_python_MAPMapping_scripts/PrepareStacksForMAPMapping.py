#%%
# this script is used to take the nrrd files that are output from the registration and downsample and smooth them in preparation for the running the "MakeTheMapMap.py" script. 
# it will output a "SmoothedTiffs" folder in the same directory as the folder containing the '.nrrd' files. you will need to separate these .tiff files into folder for the different groups and then run the the "MakeTheMapMap.py" script

import numpy as np
import nrrd
import os
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import zoom, gaussian_filter
from skimage.io import imsave
from natsort import natsorted
from PIL import Image
from tqdm.notebook import tqdm

dir = r'/media/BigBoy/ciqle/LeicaStellaris/20250225_HuCGCAMP_pERK_tERK_Ethanol_5DF/reformatted'
out_dir = dir.replace(os.path.basename(os.path.normpath(dir)), "SmoothedTiffs")
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

os.chdir(dir)

files = natsorted(glob.glob('*.nrrd'))
print(files)


for filename in tqdm(files):
    os.chdir(dir)
    IM, header = nrrd.read(filename)

    # downsample
    final_size = np.array((300, 679, 80), dtype=float)
    zoom_fact = final_size/IM.shape
    IM_resize = zoom(IM, zoom_fact)

    # gaussiuan blurr
    IM_resize = gaussian_filter(IM_resize, (2,2,2))

    if np.max(IM_resize) > 10000: # i dont remember where this number comes from. might be a problem in the future if imaging conditions change

        IM_resize = IM_resize/65535*255

    IM_resize = IM_resize.astype('uint8')

    os.chdir(out_dir)
    out_name = filename.replace('.nrrd', '_smoothed.tif')
    imsave(out_name, np.moveaxis(IM_resize, (0,1,2), (2,1,0)))


