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



HEADERS = OPS_ISEM_kernels.h OPS_ISEM_eddy_functions.h OPS_ISEM_constants.h TBL_data.h

OTHER_FILES = 
OPS_GENERATED = OPS_ISEM_ops.cpp 
OPS_FILES = OPS_ISEM.cpp 


APP=OPS_ISEM
MAIN_SRC=OPS_ISEM

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app_legacy
