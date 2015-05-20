// this simple script asks you to point to an input folder. This must contain the .nrrd warped/registered stacks output from the CMTK. It will ask you for an output folder, where the downsampled and smoothed Tiff stacks will be output. These are then processed in Matlab to creat MAP-Maps. 

max = 255; // set the max value of the images here - either 65535 for uint16, or 4096 if uint12 mapped

setBatchMode(true);


source_dir = getDirectory("Source Directory");
target_dir = getDirectory("Target Directory");

 list = getFileList(source_dir);

  for (i=0; i<list.length; i++) {
         
        open(source_dir + "/" + list[i]);        
	run("Size...", "width=300 depth=80 constrain average interpolation=Bilinear");
        run("Gaussian Blur...", "sigma=2 stack");
        setMinAndMax(0, max);
	run("8-bit");
	saveAs("tiff", target_dir + "/" + list[i] + "GauSmooth.tiff");
	close();
	showProgress(i, list.length);
	run("Collect Garbage");
  }
