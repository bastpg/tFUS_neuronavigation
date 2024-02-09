# tFUS_neuronavigation

**Objective:** We develop and disseminate a model-based navigation (MBN) tool for acoustic dose delivery in the presence of skull aberrations that is easy to use by non-specialists.   

**Please cite our paper at:**

![tFUS_gui_figures_v15](https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/46a368f3-d179-4d93-ac97-166008db5bcd)


# Methods

## Dependencies
This program requires the following depencies:
- MATLAB (version???)

## Create a Dataset
This program requires the following for each subject: 
- T1-weighted MRI Volume

## Pre-processing

```MATLAB
x = code?
```

## Pre-Calculation Walkthrough
Start the GUI by running "precalculations_GUI_Fin.mlapp" in the MATLAB Command Window

<img width="920" alt="Screen Shot 2024-02-09 at 9 35 22 AM" src="https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/f148ea3f-36ce-40d8-bae4-f0b828fda13c">

Follow the steps below,
1. Load Data
   - Click "Browse head mask" to load hmask.mat
   - Click "Browse porosity" to load poro.mat
   - Click "Browse ASEG" to aseg.mat
2. Generate Mesh
   - Click "Generate Head Mesh"
     - Modify edge length based on desire preferences (Preset of 8 mm)
   - Click "Remove invalid faces"
     - Draw a line above the eyes/ears to remove those faces from the mesh
   - Make sure normals point away from the head
     - Click "Show Normals" and ensure all red arrows point outward
     - Additionally click "Invert Normas" and ensure all red arrows point inside the mesh
3. Choose Transducer Model
   - Click a preset transducer **or**
   - Define Transducer Parameters
     - Focal Distance (mm)
     - Aperture Diameter (mm)
     - Distance to Scalp (mm)
     - Frequency (KHz)
4. Run Acoustic Intensity Calculations
   -  Pick the "Number of parallel processing units" (Preset of 10)
   -  Type the "Simulation folder name," with the full path preferred
   -  Click "Run"
  
   
     

## Pre-computation of Acoustic Beams

## Real-time tFUS Neuronavigation


# References

