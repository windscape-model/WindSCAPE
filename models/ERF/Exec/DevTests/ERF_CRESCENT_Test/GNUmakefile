# ============================================================
# ERF-CRESCENT Coupled Test Case
# DOE-funded air-sea interaction modeling (LLNL)
# ============================================================

# AMReX
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE       = FALSE
TINY_PROFILE  = FALSE
COMM_PROFILE  = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE   = FALSE
USE_GPROF     = FALSE

# Performance
USE_MPI  = TRUE
USE_OMP  = FALSE

USE_CUDA = FALSE
USE_HIP  = FALSE
USE_SYCL = FALSE

# FFT support (required for CRESCENT)
USE_FFT = TRUE

# Debugging
DEBUG = FALSE

TEST = TRUE
USE_ASSERTION = TRUE

# ===== CRESCENT COUPLING ENABLED =====
USE_CRESCENT_COUPLING = TRUE
CRESCENT_HOME = ../../../../CRESCENT

# FFTW location (required for CRESCENT coupling)
FFTW_HOME = /Users/li100/local/Cellar/fftw/3.3.10_2

# Optional: NetCDF support
# USE_NETCDF = TRUE

# GNU Make
Bpack := ./Make.package
Blocs := .

ERF_HOME := ../../..
ERF_PROBLEM_DIR = $(ERF_HOME)/Exec/DevTests/ERF_CRESCENT_Test
include $(ERF_HOME)/Exec/Make.ERF
