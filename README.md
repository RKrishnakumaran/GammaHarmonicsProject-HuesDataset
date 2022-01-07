# GammaHarmonicsProject-HuesDataset
 Repository of codes and data for replicating results of the Gamma Harmonics paper
############################################

Overview of Data:
Data collected from Microelectrode arrays implanted in V1 of macaque subjects M1 and M2 have been extracted and segmented into spikes and LFP signals and are made available within the follwoing directories:
* _SForiAchro_ directory - Microelectrode array data from M1 during fullscreen achromatic grating presentation.
* _Color_ directory - Microelectrode array data during fullscreen plain hue presentation from M1, and both fullscreen plain hue and achromatic grating presentation from M2. 

############################################

Overview of Codes: The following scripts are to be run in order to replicate the analysis and data Figures in the paper.

* _prepareData.m_ - performs Multitaper and phase analyses on LFP data for each subject and saves the output for use in Figure codes.

 <!> Requires **CHRONUX toolbox** (http://chronux.org/; tested on _version 2.00_) for Multitaper analysis 

 <!> Requires **Signal Processing Toolkit** for signal filtering and phase analysis. 

 <!> Requires **Circstat** toolbox by _Philipp Berens_ (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
  
* _Fig2.m_ - performs analysis of peak frequencies of gamma and subsequent bump and generates Fig 2

* _Fig5.m_ - performs statistical analysis of phase difference and generates Fig 5

############################################
