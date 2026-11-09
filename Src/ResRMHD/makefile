# **********************************************************************
#  Makefile for Resistive Relativitic MHD
#
#  Note: VPATH includes also the RMHD/ and RHD/ directories where
#        some shared filed are needed.
# **********************************************************************

VPATH        +=  $(SRC_HO)/ResRMHD/
VPATH        +=  $(SRC)/ResRMHD/:$(SRC)/RMHD:$(SRC)/RHD
#INCLUDE_DIRS += -I$(SRC)/ResRMHD -I$(SRC)/RHD
#INCLUDE_DIRS += -I$(SRC)/ResRMHD 
OBJ     +=  imex_source.o
OBJ     +=  mappers.o
OBJ     +=  res_eta.o rhs.o 
OBJ     +=  set_solver.o tvdlf.o mhll.o mhllc.o
HEADERS +=  resrmhd.hpp vel_limiter.hpp

# -- Files shared with other modules --

OBJ     +=  rhd_energy_solve.o rhd_pressure_fix.o \
            rmhd_energy_solve.o rhd_entropy_solve.o

# rhd_entropy_solve.o
