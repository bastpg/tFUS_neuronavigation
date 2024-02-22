# tFUS_neuronavigation

**Objective:** We develop and disseminate a model-based navigation (MBN) tool for acoustic dose delivery in the presence of skull aberrations that is easy to use by non-specialists.   

**Please cite our paper at:**

![tFUS_gui_figures_v15](https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/46a368f3-d179-4d93-ac97-166008db5bcd)


# Methods

## Dependencies
This program requires the following dependencies: MATLAB (ver. R2023a), iso2mesh Library (included in download, GitHub found [here](https://github.com/fangq/iso2mesh)), FreeSurfer, and a Linux Environment. For pre-calculation, a minimumm of 8-core CPU and 8 GB of RAM is required (≥20-core CPU and 32 GB of RAM preferred). Whereas the pre-calculation GUI should be run on a large, powerful computer, the planning GUI is best run on a local machine to avoid display lag. We recommend transferring the entire solution dataset on the local machine for a smooth viewing experience.

## T1 Pre-processing: preprocess_T1.m
This program requires a T1-weighted MRI volume for each subject. This script also uses the Freesurfer routine, which to be used within Matlab needs to be sourced on the terminal used to launch Matlab, e.g.: source /usr/local/freesurfer/nmr-dev-env (replace with your Freesurfer installation path).

Next, ensure the following lines lead to your specific tFUS library, iso2mesh library, and desired T1.

```MATLAB
close all;
clear;

addpath( '/autofs/space/guerin/USneuromod/MBN_GUI/Library' );
addpath( '/autofs/space/guerin/USneuromod/MBN_GUI/iso2mesh/' );

t1filepath = 'T1.nii.gz';  % the input T1 image (can be any format including *.mgz, *.nii, *.nii.gz etc...)
nthreads = 10;  % number of threads for SAMSEG
```

Next, run SAMSEG (~10 minutes on 10 threads), this creates the volumes seg.mgz and input/t1w/r001.mgz used hereafter. 

```MATLAB
eval( sprintf( '!samseg --t1w %s --o ./ --threads %d --refmode t1w' , t1filepath , nthreads ) );
!mri_convert input/t1w/r001.mgz T1.nii  
!mri_convert seg.mgz seg.nii
```

### Porosity Option 1: Assume uniform bone porosity (simplest, least accurate)
Assume uniform bone porosity within the SAMSEG skull mask. Ignores the pores, assumes 100% bone density everywehere in the skull.

``` MATLAB
header = load_nifti( 'seg.nii' );
smask = double( header.vol == 165);  % create skull mask from SAMSEG segmentation output
poro = 1 - double(smask);  % all skull is 100% bone
header.vol = poro;
save_nifti( header , 'poro.nii' );
```

### Porosity Option 2: Niftiweb pCT tool
Go to the  [Niftiweb pCT tool](http://niftyweb.cs.ucl.ac.uk/program.php?p=PCT), and upload the transformed T1 from SAMSEG input/t1w/r001.mgz converted to *.nii format as shown below.

```MATLAB
header = load_nifti( 'pct_0001.nii' );  % load Niftiweb pCT result
poro = 1 - header.vol / 1000;  % porosity calculation based on Aubry et al. "Experimental demonstration of noninvasive transskull adaptive focusing based on prior computed tomography scans." The Journal of the Acoustical Society of America 113.1 (2003): 84-93.
poro( find(poro>1) ) = 1;  % ensure porosity is between 0 and 1
poro( find(poro<0) ) = 0;
header.vol = poro;
save_nifti( header , 'poro.nii' )
```

### Porosity Option 3: mri-to-ct deep learning tool
mri-to-ct details can be found [here](https://github.com/MatDagommer/mri-to-ct). The first step is to compute a liberal mask of the skull as shown below (input needed to mri-to-ct)

```MATLAB
header = load_mgh( 'seg.nii' );
smask = double( header.vol == 165);  % create skull mask from SAMSEG segmentation output
smask = imclose( smask , strel('sphere',5) );  % close the skull mask
smask = thickenbinvol( smask , 3 );  % make sure the mask covers a bit more than the actual skull (liberal mask)
header.vol = smask;
save_nifti( header , 'skull_mask.nii' );  % save skull mask as NIFTI to be used as input to mri-to-ct
```

After creating the porosity nifti (poro.nii), run the following lines to create the porosity, head mask, and aseg files.

```MATLAB
header = load_mgh( 'seg.nii' );
hmask = header.vol > 0;
hmask = imclose( hmask , strel('sphere',5) );  % close the mask, fill holes and smooth
for ii = 1 : size(hmask,3)
    hmask(:,:,ii) = imfill( hmask(:,:,ii) , 'holes' );
end
hmask = smooth3( double(hmask) , 'box' , 3 ) > 0.5;
header.vol = hmask;
save_nifti( header , 'hmask.nii' )

A = load_nifti( 'poro.nii' );
B = load_nifti( 'hmask.nii' );
C = load_nifti( 'aseg.nii' );

poro = reslice_( A.vol );
hmask = reslice_( B.vol );
aseg = reslice_( C.vol );

save  poro_resliced.mat  poro;
save  aseg_resliced.mat  aseg;
save  hmask_resliced.mat  hmask;
```




## Pre-Calculation: precalculations_GUI_Fin.mlapp
Start the GUI by running "precalculations_GUI_Fin" in the MATLAB Command Window

<img width="920" alt="Screen Shot 2024-02-09 at 9 35 22 AM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/f148ea3f-36ce-40d8-bae4-f0b828fda13c">

### Step 1: Load Data
First, load the data. The indicators next to each button should turn green after the file is loaded. Then, use the display to check for any inconsitences or mistakes with the headmask, porosity, and aseg files.

<img width="763" alt="Screen Shot 2024-02-15 at 6 56 24 AM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/8a4d9566-88ab-4075-8946-18d0a7fb5660">


### Step 2: Generate Mesh
This step generates the mesh. The edge length is the average length of triangle discretization of the scalp mesh. In other words, a small edge length (4 mm) yields dense discrete triangles while a large edge length (16 mm) yields a coarse mesh. We have found 8 mm balances the accuracy of a dense mesh with computational efficiency of a course mesh. After generating the mesh, remove the faces below and including the ears and eyes. This should ideally yield a mesh with 2000-3000 faces. Lastly, check that all normals are pointing away from the mesh. 


<img width="1405" alt="Screen Shot 2024-02-15 at 7 05 35 AM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/ea73e972-0f57-4b3d-a395-49e465f9fb94">


### Step 3: Choose Transducer Model
In this step, either select a pre-existing transducer model or design custom parameters.

<img width="1020" alt="Screen Shot 2024-02-15 at 7 41 00 AM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/a745c355-e311-43e5-923e-5dfe423f8e96">


### Step 4: Run Acoustic Intensity Calculations
In this step, pick the number of parallel processing units and indicate the file path. For parallel processing units, we recommend 20 parallel units which gives a pre-computation time of ~0.5 hours for 2500 faces. Computational time will increase with the density of the scalp mesh. Lastly, specify the absolute path for the output folder (using / to separate directories). Then, click "Run."

<img width="757" alt="Screen Shot 2024-02-15 at 7 52 24 AM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/255aa677-6050-48db-9e42-b3a910e23b4b">


### Step 5: Check the Data Output
The output folder should include four .mat files: app_data, mSOUND_inputs, SOL_Beams, SOL_ScalpMaps. app_data.mat and mSOUND_inputs.mat are both input files into the pre-calculation. SOL_Beams.mat is the largest file (a few GB depending on the density of the mesh) and stores the beam calculations. SOL_ScalpMaps.mat includes the scalp maps generated from the beams.
     

## Neuronavigation Planning GUI: planning_GUI_Fin.mlapp
Whereas the pre-calculation GUI should be run on a large, powerful computer, the planning GUI is best run on a local machine to avoid display lag. We recommend transferring the entire solution dataset on the local machine for a smooth viewing experience. 

First, load the folder with the pre-calculation results. You must load the entire folder, rather than individual .mat files. On the bottom of the GUI in red text, the progress of loading the pre-calculation results will be indicated. Then, select the desired nuclei and click "Display scalp map." 

<img width="1063" alt="Screen Shot 2024-02-18 at 12 11 48 PM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/864c427b-7094-43c4-a1db-fa3180db5160">

The main display window is interactive. The scalp map can me moved either using pre-determined options on the left under "NAVIGATE" or using traditional MATLAB figure controls on the top right of the display. Nuclei visibility can be modified in the bottom left with the light blue nuclei the current target. The head mesh can also be displayed along with only showing a single side of the scalp map. 

<img width="1406" alt="Screen Shot 2024-02-18 at 12 19 34 PM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/a8e198d5-d003-4d38-9582-f0fd9d35a694">

Additionally, you can click on a scalp map position and the display upbates the transducer, beam, and intensity bar graph at the new location. MATLAB tools in the upper right corner must be inactivated in order to change the transducer location. For easier positioning of the coil, it is best to click on the outside of the scalp map.

<img width="1394" alt="Screen Shot 2024-02-18 at 12 28 38 PM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/f8ce2003-f64b-40d2-ba5d-9a18911153e7">

The power distribution and scalp map figures can be saved using the "Save figures" button. This saves two figures in a single Matlab file, which can then be loaded using the load command. The .mat file saved will be located in the originally loaded folder with the following naming convention: FIG_*NUCLEI*_transd*LOCATIONNUMBER*.mat. For example, a saved figure for the left thalamus at trandsducer location #1505 will be saved as: FIG_LeftThalamus_transd1505.mat. Unlike in the GUI, the transducer location cannot be changed.

<img width="1163" alt="Screen Shot 2024-02-18 at 12 38 07 PM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/a7fdcf5a-07d4-4cf2-9833-243818e56337">

The visibility of objects in the scalp map figure can be changed by opening the plot browser in View > Plot Browser:

<img width="1443" alt="Screen Shot 2024-02-18 at 12 41 48 PM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/8fe1b0f8-890c-479f-b7fa-faa5697bcbf4">

