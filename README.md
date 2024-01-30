# tFUS_neuronavigation

**Background:** Transcranial focused ultrasound (tFUS) neuromodulation has shown promise in animals but is hard to translate to humans because of the thicker skull that absorbs and scatters ultrasound waves. 

**Objective:** We develop and disseminate a model-based navigation (MBN) tool for acoustic dose delivery in the presence of skull aberrations that is easy to use by non-specialists.   

**Methods:** We pre-compute acoustic beams for thousands of virtual transducer locations on the scalp of the subject under study. We use the hybrid angular spectrum solver mSOUND which runs in ~4 seconds per solve per CPU, yielding pre computation times under one hour for scalp meshes with up to 4,000 faces and a parallelization factor of 5. We combine this pre-computed set of solutions with optical tracking of the transducer, allowing real-time display of the tFUS beam as the operator freely moves the transducer. We assess the impact of MBN versus line-of-sight targeting (LOST) positioning in 13 test subject simulations.

**Results:** Our navigation tool has a display refresh rate of ~2 Hz. In our simulations, MBN increased the acoustic dose in the thalamus and amygdala of the test subjects by 22-137% compared to LOST and avoided complete target misses that affected 10-20% of LOST cases. MBN yielded a lower variability of the deposited dose across subjects than LOST.

**Conclusions:** MBN may yield greater and more consistent (less variable) ultrasound dose deposition than LOST, and therefore has the potential to improve the significance and variability of outcome measures of tFUS neuromodulation. 

**Please cite our paper at:**

![tFUS_gui_figures_v15](https://github.com/parkerkotlarz/tFUS_neuronavigation/assets/157265957/46a368f3-d179-4d93-ac97-166008db5bcd)


# Methods

## Dependencies

## Create a Dataset


## Pre-processing

```MATLAB
x = code?
```

## GUI Walkthrough

## Pre-computation of Acoustic Beams

## Real-time tFUS Neuronavigation


# References

