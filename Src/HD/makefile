# **********************************************************************
# Makefile for the HydroDynamics (HD) Module
#
# Note that we add Src/MHD to VPATH since HD and MHD
# share some files (see below)
# **********************************************************************

VPATH        +=  $(SRC)/HD:$(SRC)/MHD
INCLUDE_DIRS += -I$(SRC)/HD
OBJ          += mappers.o hll.o hllc.o \
                rhs.o roe.o set_solver.o tvdlf.o two_shock.o
HEADERS      += hd.hpp
