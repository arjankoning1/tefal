# TEFAL

TEFAL is a software package for the processing of the output files of the TALYS nuclear reaction code, and data from other sources, into an ENDF-6 formatted nuclear data library.

## Documentation and reference

A description of the code and its options can be found in the [TEFAL Tutorial (pdf)](https://github.com/arjankoning1/tefal/blob/main/doc/tefal.pdf).

The reference to be used for TEFAL is:

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155, 1 (2019).

## Installation

### Prerequisites

The following are the prerequisites for compiling TEFAL:

- GNU make
- a recent Fortran compiler, such as GNU Fortran (gfortran)
- a successful installation of the TALYS nuclear model code
- git, only when TEFAL is downloaded using `git clone`

### Downloads

TEFAL can be downloaded in one of the following ways.

#### 1. Frozen version TEFAL-2.2 (December 2025)

The frozen TEFAL-2.2 distribution is available from the [TALYS page](https://nds.iaea.org/talys/). It can be retrieved by clicking on the download link or with

```bash
curl -LO https://nds.iaea.org/talys/codes/tefal.tar
tar zxf tefal.tar
```

This version is fixed and will not change.

#### 2. Latest beta version without git

Users who do not have git can download a snapshot of the current `main` branch directly from GitHub:

```bash
curl -L \
  -o tefal-main.tar.gz \
  https://github.com/arjankoning1/tefal/archive/refs/heads/main.tar.gz

tar zxf tefal-main.tar.gz
mv tefal-main tefal
```

This produces the same `tefal/` directory structure as the git version, but without the git history.

The downloaded snapshot contains the latest version of the `main` branch at the time of download. To obtain a newer version later, download the snapshot again.

#### 3. Latest beta version using git

Users with git can clone the repository with

```bash
git clone https://github.com/arjankoning1/tefal.git
```

The advantage of this method is that the local TEFAL installation can subsequently be updated with

```bash
cd tefal
git pull --ff-only
```

### Installation instructions

#### 1. Frozen version TEFAL-2.2

For the frozen tar distribution:

```bash
cd tefal
./install_tefal.bash
```

An alternative option is:

```bash
cd tefal/source
make
```

The above will invoke the default compiler `gfortran`.

#### 2. Latest beta version

The installation procedure is identical whether the latest beta version was obtained as a GitHub tar snapshot or using `git clone`.

From the `tefal/` directory, run

```bash
./install_tefal.bash
```

which automatically executes the `Makefile` in `tefal/source`. At the end, `install_tefal.bash` will print the recommended shell configuration.

An alternative option is:

```bash
cd tefal/source
make
```

The executable is installed as

```text
tefal/bin/tefal
```

For the latest beta version, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, the Makefile uses:

```text
-w -O3 -ffp-contract=off
```

For other compilers, no default compiler flags are imposed.

The compiler and compilation options can be passed to the Makefile through `install_tefal.bash`. For example:

```bash
# GNU Fortran
./install_tefal.bash FC=gfortran FFLAGS="-O3 -ffp-contract=off"

# Intel Fortran
./install_tefal.bash FC=ifx FFLAGS="-O3"
```

Set `TEFAL_DIR` to the TEFAL installation directory. This variable is required unless the fallback path in `source/machine.f90` has been set manually. For example:

```bash
export TEFAL_DIR="/Users/koning/tefal"
```

If you want to run `tefal` from anywhere, add the TEFAL `bin` directory to `PATH`:

```bash
export PATH="$TEFAL_DIR/bin:$PATH"
```

These lines can be added to your shell configuration file, for example `~/.zshrc` or `~/.profile`.

If setting `TEFAL_DIR` is not possible on a particular system, edit `code_dir` in `source/machine.f90` and rebuild TEFAL.

## The TEFAL package

The `tefal/` directory contains the following directories and files:

- `README.md` this README file
- `LICENSE` the License file
- `install_tefal.bash` installation script
- `source/` the Fortran source code of TEFAL and the Makefile
- `bin/` the `tefal` executable after successful installation
- `misc/` text files and energy grids used by TEFAL
- `doc/` the tutorial in PDF format
- `samples/` the input and output files of the sample cases, and the `verify` script used to run the sample cases

In total, about 600 MB of free disk space is required to install TEFAL.

## Sample cases

The sample cases provide examples of the use of TEFAL and can be used to verify a successful installation. The `samples/` directory contains various sample cases with a subdirectory `org/` containing the reference results and a subdirectory `new/` for results produced by the user.

The TEFAL sample cases assume that TALYS is installed in a sibling directory:

```text
parent_directory/
├── talys/
└── tefal/
```

To run the sample cases:

```bash
cd samples
./verify
```

From the top-level TEFAL directory, the same test can be started with:

```bash
make -C source check
```

`make check` sets `TEFAL_DIR` and the sibling `TALYS_DIR` automatically for the test.

TEFAL normally operates on the result of a TALYS calculation performed with the keyword:

```text
endf y
```

For example, with `talys.inp` and `tefal.inp`:

```bash
talys < talys.inp > talys.out
tefal < tefal.inp > tefal.out
```

assuming that both `talys/bin` and `tefal/bin` are available through `PATH`.

## License and Copyright

This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
