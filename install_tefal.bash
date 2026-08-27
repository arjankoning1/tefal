#!/usr/bin/env bash

set -euo pipefail

# Determine the TEFAL installation directory, independently of where
# the script is called from.

tefal_dir=$(cd "$(dirname "$0")" && pwd)
source_dir="$tefal_dir/source"

# Verify that the expected directories and build files exist.

if [[ ! -d "$source_dir" ]]; then
  echo "TEFAL installation error: source directory not found:" >&2
  echo "  $source_dir" >&2
  exit 1
fi

if [[ ! -f "$source_dir/Makefile" ]]; then
  echo "TEFAL installation error: Makefile not found:" >&2
  echo "  $source_dir/Makefile" >&2
  exit 1
fi

misc_file="$tefal_dir/misc/endf_n.txt"

if [[ ! -f "$misc_file" ]]; then
  echo "TEFAL installation error: misc database missing or incomplete:" >&2
  echo "  $misc_file" >&2
  exit 1
fi

echo
echo "Installing TEFAL"
echo "Installation directory: $tefal_dir"
echo

# Pass all command-line arguments directly to make. This permits, e.g.:
#
# ./install_tefal.bash FC=ifx FFLAGS="-O3"
# ./install_tefal.bash FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"

make -C "$source_dir" clean
make -C "$source_dir" all "$@"

tefal_exe="$tefal_dir/bin/tefal"

if [[ ! -x "$tefal_exe" ]]; then
  echo "TEFAL installation error: executable not created:" >&2
  echo "  $tefal_exe" >&2
  exit 1
fi

echo
echo "TEFAL executable:"
echo "  $tefal_exe"
echo
echo "If not already done, add the following lines to your shell configuration:"
echo
echo "  export TEFAL_DIR=\"$tefal_dir\""
echo "  export PATH=\"\$TEFAL_DIR/bin:\$PATH\""
echo
echo "Alternatively, edit code_dir in source/machine.f90 and rebuild TEFAL."
echo
