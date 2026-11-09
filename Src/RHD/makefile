# **********************************************************************
# Makefile for the Relativistic HydroDynamics (RHD) Module
#
# Note that we add Src/RMHD to VPATH since RHD and RMHD
# module share some files (see below)
# **********************************************************************

VPATH        +=  $(SRC)/RHD:$(SRC)/RMHD
INCLUDE_DIRS += -I$(SRC)/RHD
OBJ          += hll.o hllc.o mappers.o set_solver.o \
                rhd_energy_solve.o rhd_entropy_solve.o rhd_pressure_fix.o \
                rhs.o tvdlf.o two_shock.o
HEADERS      += rhd.hpp vel_limiter.hpp

