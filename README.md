
# TEFAL
TEFAL is a software package for the processing of the output files of the TALYS nuclear reaction code, and data from other sources, into an ENDF-6 formatted nuclear data library.

## Documentation and reference
A description of the code and its options can be found in the [TEFAL Tutorial (pdf)](https://github.com/arjankoning1/tefal/blob/main/doc/tefal.pdf).
The reference to be used for TEFAL is

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155,1 (2019).

## Installation

### Prerequisites:

The following are the prerequisites for compiling TEFAL:
  - git (only if the package is downloaded via Github)
  - a recent Fortran compiler such as gcc (gfortran)
  - a successful installation of the TALYS nuclear model code

### Downloads:

To download TEFAL, you can use one of the following options:
#### 1. Download the entire tar file (frozen version):
```
https://nds.iaea.org/talys/tefal.tar
tar zxf tefal.tar
```

#### 2. Using git (latest beta version):
```
git clone https://github.com/arjankoning1/tefal.git
```
The TEFAL sample cases do not fall under the git repository. For that you need to download:
```
https://nds.iaea.org/talys/samples/tefal_samples.tar
tar zxf tefal_samples.tar
```
The resulting *samples/* directory should be moved inside the *tefal/* directory.

### Installation instructions:

To install TEFAL, you can use one of the following options:
#### 1. Using make:
```
cd tefal/source
make
```
#### 2. Using the install_tefal.bash script:
```
cd tefal
install_tefal.bash tefal
```

The above will produce a *tefal* executable in the *tefal/bin* directory. 
The compiler and its flags can be set in either the *source/Makefile* or in *code_build*.

## The TEFAL package

The *tefal/* directory contains the following directories and files:

+ `README.md` this README file
+ `LICENSE` the License file
+ `install_tasman.bash`, `code_build.bash` and `path_change.bash` installation scripts
+ `source/` the Fortran source code of TEFAL and the Makefile
+ `bin/` the executable after successful installation
+ `misc/` text files and energy grids to be used by TEFAL
+ `doc/` the tutorial in pdf format
+ `samples/` the input and output files of the sample cases, and the *verify* script to run the sample cases

In total, you will need about 600 Mb of free disk space to install TEFAL.

## Sample cases

The sample cases serve to provide examples of the use of TEFAL and to verify a successful installation. The *samples/* directory contains various sample cases with a subdirectory *org/* with our results and a subdirectory *new/* with the results produced by the user. The entire sample set will take about 1 hour.
```
cd samples
./verify
```

TEFAL will only work after a TALYS calculation with the keyword **endf y** in the TALYS input file.
You may create your own input files, e.g. *talys.inp* and *tefal.inp*, after which TEFAL works as follows:

```
talys < talys.inp > talys.out
tefal < tefal.inp > tefal.out
```
assuming that *tefal* is linked to the *tefal/bin/tefal* executable and *talys* is linked to the *talys/bin/talys* executable.

## License and Copyright
This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
