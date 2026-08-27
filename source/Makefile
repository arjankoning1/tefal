# TEFAL Makefile
#
# Examples:
#
#   make
#   make FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"
#   make FC=ifx FFLAGS="-O3"
#   make clean
#   make check

# Compiler-specific flags can be supplied by the user.
ifeq ($(origin FC), default)
  FC = gfortran
endif

ifeq ($(origin FFLAGS), undefined)
  ifneq ($(findstring gfortran,$(FC)),)
    FFLAGS = -w -O3 -ffp-contract=off
  else
    FFLAGS =
  endif
endif

# Installation directories.

TEFAL_DIR := $(abspath ..)
BIN_DIR   := $(TEFAL_DIR)/bin
TARGET    := $(BIN_DIR)/tefal

TALYS_DIR_CHECK := $(abspath $(TEFAL_DIR)/../talys)

# TEFAL source files.
#
# Free-form source files must be named so that alphabetical order places
# module definitions, such as A0_tefal_mod.f90, before their users.

SOURCES := $(sort $(wildcard *.f) $(wildcard *.f90))

.PHONY: all install clean check

all install: $(TARGET)

$(TARGET): $(SOURCES)
	@echo "Compiler: $(FC)"
	@echo "Flags:    $(FFLAGS)"
	@mkdir -p "$(BIN_DIR)"
	$(FC) $(FFLAGS) -o "$(TARGET)" $(SOURCES)
	@rm -f *.o *.mod
	@echo "Installed: $(TARGET)"

clean:
	rm -f *.o *.mod tefal "$(TARGET)"

check: all
	cd "$(TEFAL_DIR)/samples" && \
	  TEFAL_DIR="$(TEFAL_DIR)" \
	  TALYS_DIR="$(TALYS_DIR_CHECK)" \
	  ./verify
