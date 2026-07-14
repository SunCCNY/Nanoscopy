# SMLM_Lib

SMLM_Lib is a Python library in SMLM developed in the Nanoscopy Laboratory at the City College of New York. 

SMLM_Lib can be broadly applied to modeling, theoretical analysis, simulation, image reconstruction, performance benchmarking and evaluation in SMLM. 

SMLM_Lib is used in our projects in SMLM. 

SMLM_Lib is built upon the universal model of frames [1, 5], the Markov chain model of fluorescence molecule photoactivation process [2, 4], and quality metrics [3, 5]. 

## Functions

SMLM_Lib includes the following data classes:
- Camera: all types [1]
- Poisson noise: uniform [1] and non-uniform
- Gaussian noises: uniform [1] and pixel-variant for sCMOS cameras
- 2D PSF: Gaussian [1], Airy [1]
- 3D PSF: Astigmatic [1]
- Emitter photoactivation processes: Markov chain model [2, 8]

and functions for above data classes: 
- Frame generator [1] 
- Movie generator [2, 3, 8]
- SNR: uniform noise [4, 6], non-uniform noise 
- Quality metrics: root mean square minimum distance (RMSMD) [2], RMSMD with partition (RMSMD-P) [2], and root mean square error (RMSE) with partition (RMSE-P) [2]
- Fisher information matrix [1, 5]
- unbiased Gaussian information-achieving estimator for a frame (UGIA-F) [2]
- EM algorithm 


## Acknowledgements

This work was partly supported by the NSF under Grant Number CCF-2313072 and the Army Research Office under Grant Number W911NF-23-1-0189.

The views and conclusions contained in this document are those of the authors and should not be interpreted as representing the official policies, either expressed or implied, of the Army Research Office or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for Government purposes notwithstanding any copyright notation.

## References

[8] Y. Sun, "Markov chain models of emitter activations in single molecule localization microscopy," Optics Express, 32(19), 33779-33794(2024).

[7] Y. Sun, "Partition of estimated locations: an approach to accurate quality metrics for stochastic optical localization nanoscopy," JOSA A, 39(12), 2307-2315(2022).

[6] M. Sun and Y. Sun, "Information sufficient segmentation and signal-to-noise ratio for 3D astigmatism stochastic optical localization nanoscopy," Electr. Letters, 58(2), 58-60(2022). 

[5] Y. Sun, and Y. Guan, "Effect of unknown emitter intensities on localization accuracy in stochastic optical localization nanoscopy using single frames," JOSA A, 38(12), 1830-1840(2021). 

[4] Y. Sun, "Information sufficient segmentation and signal-to-noise ratio in stochastic optical localization nanoscopy," Optics Letters, 45(21), 6102-6105(2020). 

[3] Y. Sun, "Spatiotemporal resolution as an information theoretical property of stochastic optical localization nanoscopy," 2020 Quantitative BioImaging Conf. (QBI2020), Oxford, UK, Jan. 6-9, 2020. 

[2] Y. Sun, "Root mean square minimum distance as a quality metric for stochastic optical localization nanoscopy images," Sci. Reports, 8(1), 17211(2018). 

[1] Y. Sun, "Localization precision of stochastic optical localization nanoscopy using single frames," J. Biomed. Optics, 18(11), 111418-14(2013). 

## Contact
Yi Sun, Electrical Engineering Department, Nanoscopy Laboratory, The City College of City University of New York, New York, NY 10031, USA. E-mail: ysun@ccny.cuny.edu

