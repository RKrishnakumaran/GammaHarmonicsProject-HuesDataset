# GammaHarmonicsProject-HuesDataset
 Repository of codes and data for replicating results of the Gamma Harmonics paper
############################################

Overview of Data:
Data collected from Microelectrode arrays implanted in V1 of macaque subjects M1 and M2 have been extracted and segmented into spikes and LFP signals and are made available within the follwoing directories:
* _SForiAchro_ directory - Microelectrode array data from M1 during fullscreen achromatic grating presentation.
* _Color_ directory - Microelectrode array data during fullscreen plain hue presentation from M1, and both fullscreen plain hue and achromatic grating presentation from M2. 

############################################

Overview of Codes: The following scripts, run in order, would replicate the gamma-harmonic analyses and Figures 2 and 5 in the paper from the LFP data.

* _prepareData_r2019b.m (for MATLAB R2019b release)_ or _prepareData.m (for future releases, in case findpeaks returns a structure)_ - performs Multitaper and phase analyses on LFP data for each subject and saves the output for use in Figure codes.

 <!> Requires **CHRONUX toolbox** (http://chronux.org/; tested on _version 2.00_) for Multitaper analysis 

 <!> Requires **Signal Processing Toolkit** for signal filtering and phase analysis. 

 <!> Requires **Circstat** toolbox by _Philipp Berens_ for directional statistics. (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
  
* _Fig2.m_ - performs analysis of peak frequencies of gamma and subsequent bump and generates Fig 2

* _Fig5.m_ - performs statistical analysis of phase difference and generates Fig 5

############################################
