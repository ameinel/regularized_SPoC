# regularized_SPoC
This repository contains the Matlab implementation of different regularization strategies for the SPoC algorithm:

> A. Meinel, J.S. Castaño-C., B. Blankertz, F. Lotte, M. Tangermann, "Characterizing Regularization Techniques for Spatial Filter Optimization in Oscillatory EEG Regression Problems", Springer Neuroinformatics, 2018 - [link](http://dx.doi.org/10.1007/s12021-018-9396-7)

The original publication of the SPoC algorithm without regularization can be found here:

> S. Dähne, F. C. Meinecke, S. Haufe, J. Höhne, M. Tangermann, K. R. Müller, V. V. Nikulin, "SPoC: a novel framework for relating the amplitude of neuronal oscillations to behaviorally relevant parameters", NeuroImage, 86:111-122, 2014

The original implementation can be found here: [matlab_SPoC](https://github.com/svendaehne/matlab_SPoC)

## Getting started
1. Download the following toolboxes and save them under ./external:
* [bbci_toolbox](https://github.com/bbci/bbci_public) - Matlab toolbox for BCI experimenting
* [Post-HocLabeling](https://github.com/bsdlab/Post-HocLabeling) - Framework to generate labeled data sets from arbitrary, pre-recoreded EEG files
2. Download EEG data available at https://zenodo.org/record/1065107#.WhaCX3XyvCI and place the .mat files under ./data
3. Run the script "test_regularizedSPoC.m" to get familiar with the usage in an example decoding scenario.

___

@Andreas Meinel, 06.08.2018. 
Contact: andreas.meinel@blbt.uni-freiburg.de
