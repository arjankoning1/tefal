# TEFAL

TEFAL is a software package for the processing of the output files of the TALYS nuclear reaction code, and data from other sources, into an ENDF-6 formatted nuclear data library.

## Documentation and reference

A description of the code and its options can be found in the [TEFAL Tutorial (pdf)](https://github.com/arjankoning1/tefal/blob/main/doc/tefal.pdf).

The reference to be used for TEFAL is:

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155, 1 (2019).

## Installation

### Prerequisites

The following are the prerequisites for compiling TEFAL:

- git (only if the package is downloaded via GitHub)
- GNU make
- a recent Fortran compiler, such as GNU Fortran (gfortran)
- a successful installation of the TALYS nuclear model code

### Downloads

To download TEFAL, you can use one of the following options.

#### 1. Download the entire tar file (frozen version TEFAL-2.2)

```bash
curl -LO https://nds.iaea.org/talys/tefal.tar
tar zxf tefal.tar
```

#### 2. Using git (latest beta version)

```bash
git clone https://github.com/arjankoning1/tefal.git
```

### Installation instructions

#### 1. For the tar file (frozen version TEFAL-2.2)

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

#### 2. For the git version (latest beta version)

```bash
cd tefal
./install_tefal.bash
```

which automatically executes the `Makefile` in `tefal/source`. At the end, `install_tefal.bash`
will print the recommended shell configuration.

An alternative option is:

```bash
cd tefal/source
make
```

For both installation methods, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, the Makefile uses:

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

The above will produce a `tefal` executable in the `tefal/bin` directory.

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
.../talys/
.../tefal/
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

`make check` sets `TEFAL_DIR` automatically for the test.

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
