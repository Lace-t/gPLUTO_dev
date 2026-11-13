# **********************************************************************
#  Makefile for Ideal Relativistic MHD
#
#  Note: VPATH includes also the RMHD/ and RHD/ directories where
#        some shared filed are needed.
# **********************************************************************

VPATH        +=  $(SRC)/RMHD
INCLUDE_DIRS += -I$(SRC)/RMHD -I$(SRC)/RHD
OBJ          += pcp_nr.o pcp_stage.o rmhd_consistency.o mappers.o set_solver.o hll.o hllc.o hlld.o \
                rmhd_energy_solve.o \
                rmhd_pressure_fix.o rmhd_entropy_solve.o \
                tvdlf.o rhs.o gforce.o
HEADERS      += rmhd.hpp vel_limiter.hpp

pcp_nr.o: $(SRC)/RMHD/pcp_math.hpp
