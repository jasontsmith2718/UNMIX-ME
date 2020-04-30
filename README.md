# UNMIX-ME (Unmix Multiple Emissions)
Deep Learning for Hyperspectral Fluorescence Lifetime Image (HFLI) Unmixing

--------------------------------------------------------------

This GitHub contains relevant script, data and instructions for:
1. FLI data simulation workflow (**MATLAB**)

  * [HMFLI](https://github.com/jasontsmith2718/UNMIX-ME/tree/master/dataSimulation/twoSpectraFourLifetime) (Hyperspectral Macroscopic Fluorescence Lifetime Imaging)
  * [HMFLI-RET](https://github.com/jasontsmith2718/UNMIX-ME/tree/master/dataSimulation/FRET) (HMFLI with Resonance Energy Transfer correction). **Validated for AF700/AF750 FRET pair**
  * [Other](https://github.com/jasontsmith2718/UNMIX-ME/blob/master/miscellaneous.m) - for plotting hyperspectral TPSF information and ex. generation of emission spectral data.

2. Neural network training (**python**, [_Tensorflow & Keras_]).

  * **UNMIX-ME** [CNN model](https://github.com/jasontsmith2718/UNMIX-ME/tree/master/UNMIX-ME_CNN) for HMFLI analysis

Authors: [Jason T. Smith](https://www.researchgate.net/profile/Jason_Smith96), [Marien Ochoa](https://scholar.google.com/citations?user=CiT-IycAAAAJ&hl=es)

--------------------------------------------------------------

### Relevant Work:

#### [Original preprint](https://www.biorxiv.org/content/10.1101/745216v2)
- Smith JT, Ochoa M, Intes X. "Spectral and lifetime fluorescence unmixing via deep learning." bioRxiv (2020): 745216.

#### [OSA Biophotonics Congress Proceeding & Virtual Presentation](https://www.osapublishing.org/abstract.cfm?uri=OTS-2020-SW1D.5)
- Smith JT, Ochoa M, Intes X. _"Hyperspectral Lifetime Unmixing via Deep Learning,_" in Biophotonics Congress: Biomedical Optics 2020 (Translational, Microscopy, OCT, OTS, BRAIN), OSA Technical Digest (Optical Society of America, 2020), paper SW1D.5

--------------------------------------------------------------
