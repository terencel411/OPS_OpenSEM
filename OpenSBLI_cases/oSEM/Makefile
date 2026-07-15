#
# The following environment variables should be predefined:
#
# OPS_INSTALL_PATH
# OPS_COMPILER (gnu,intel,etc)
#

USE_HDF5=1
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5

CXXFLAGS += -std=c++17


HEADERS = opensbliblock00_kernels.h defdec_data_set.h 3D_slice_tracker.hh bmp.hh

OTHER_FILES = bmp.cc
OPS_GENERATED = opensbli_ops.cpp 
OPS_FILES = opensbli.cpp 


APP=opensbli
MAIN_SRC=opensbli

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app_legacy
