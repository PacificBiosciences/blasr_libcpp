[![Build Status](https://travis-ci.org/PacificBiosciences/blasr_libcpp.svg?branch=master)](https://travis-ci.org/PacificBiosciences/blasr_libcpp)
#What is blasr_libcpp#

Blasr_libcpp is a **library** used by blasr and other executables such as samtoh5, loadPulses for analyzing PacBio sequences. This library contains three sub-directories, including pbdata, hdf and alignment.
+ pbdata  
   - contains source code for handling Pacbio sequences and can build lib ```libpbdata```.
+ hdf 
   - contains source code for handling Pacbio hdf5 files (e.g., *.pls.h5, *.rgn.h5, *.bas.h5) and builds ```libpbhdf```.
+ alignment 
   - contains source code for aligning Pacbio reads to target sequences used in blasr and builds ```libblasr```.


##Appendix: Dependencies##
+ libpbdata
   -  does **not** depend on ```libhdf5```
   -  should build without the ```pbbam``` library *for now*

- libpbhdf
   -  depends on the ```libpbdata``` and ```libhdf5``` libraries to build
   
- alignment
   -  depends on the ```libpbdata``` library to build
   -  can build either with or without the ```libpbhdf``` library
   -  can build either with or without the ```pbbam``` library
