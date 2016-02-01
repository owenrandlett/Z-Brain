# Z-Brain

This repository contains the analysis scripts and viewer function for the Z-Brain and MAP-Mapping methods from Randlett et al., 2015 ::  http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3581.html

The Z-Brain

The Z-Brain was built as a neuroanatomical reference atlas for the zebrafish neuroscience community. It was built upon a 6dpf Nacre/mitfa mutant larvae stained with anti-ERK. By registering fish stained with different anatomical markers (antibody stains, transgenes, dye fills) to this common reference brain, we have created a platform where many labels can be explored in the same reference space. We currently have 29 labels registered to the atlas. For each label we register multiple fish, and then represent that label as the mean across fish, allowing for visualization of the average positioning and staining of these neurons or features. Based on these markers, we have also segmented to brain into 294 regions.

We hope that this resource will be valuable for exploring and understanding the 6dpf zebrafish brain, and for automated and standardized annotation of new anatomical and functional data. This can be accomplished by registering new datasets into the atlas, allowing for direct comparison to the accumulated dataset. We encourage anyone interested in helping us expand and improve the atlas, either by incorporating additional labels, or by expanding or improving the regional segmentation, to please contact Owen Randlett at zebrafishbrain@gmail.com.

Using the Z-Brain to explore anatomy

You can do some simple in-browser exploration of the database at our website: http://engertlab.fas.harvard.edu/Z-Brain. 
To fully explore the Z-Brain, you will need to download the anatomy and mask databases, and run ZBrainViewer.m in Matlab. This will allow you to perform multicolour overlays of labels and regions, ‘click to find’ regions, view in different slice orientations, write slices or stacks to disk, and draw or update regional definitions. Instructions for the various functions are displayed when the program starts, and commented in the header of the ‘ZBrainViewer.m’ function.

The three files you will need are:

ZBrainViewer.m
- Contained in this github repository.

Anatomy Label Database
- http://engertlab.fas.harvard.edu/zDownloads/AnatomyLabelDatabase.hdf5
- This is an image file that contains the stacks of the labels (transgenes, antibody stains, etc). Each stack is a mean across multiple fish. It is an HDF5 file, and so can also be loaded into different viewing programs. The resultion of the data is x/y/z = 0.798/0.798/2um. 

Warning: file size is very large (~4.5gb).


Mask Database
- http://engertlab.fas.harvard.edu/zDownloads/MaskDatabase.mat
- This is the file that contains the regional definitions. It is compressed a Matlab file. If you would like the regions in an image format, please contact me.
Download these three files and put them in your Matlab path. Then run the ‘ZBrainViewer’ function.

MAP-Mapping

In order to create MAP-Maps from pERK and tERK stained data, data must first be registered to the Z-Brain reference brain stack (http://engertlab.fas.harvard.edu/zDownloads/Ref20131120pt14pl2.nrrd). To do this we use CMTK (https://www.nitrc.org/projects/cmtk/).

Once the data is registered, you can use the 'MakeTheMAPMap.m' function to create the activity map. Then the 'ZBrainAnalysisOfMAPMaps.m' function can be used to analyze MAP-Maps using the Z-Brain. 

Reference:
