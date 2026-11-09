# **********************************************************************
#  Makefile for Resistive Relativitic MHD
#
#  Note: VPATH includes also the RMHD/ and RHD/ directories where
#        some shared filed are needed.
# **********************************************************************

VPATH        +=  $(SRC)/Radiation/:$(SRC)/RMHD:$(SRC)/RHD
INCLUDE_DIRS += -I$(SRC)/Radiation
OBJ     +=  rad_hll.o rad_hllc.o rad_tvdlf.o rad_set_solver.o rad_step.o  rad_rk_step.o
OBJ     +=  rad_update_stage.o rad_rhs.o rad_substepping.o irradiation.o
HEADERS      += radiation.hpp 
#rad_tools.hpp flux_lim.hpp rad_errors.hpp rad_sourcef.hpp rad_speed.hpp newtonsolve.hpp


# -- Files shared with other modules --

# OBJ     +=  rhd_energy_solve.o rhd_pressure_fix.o \
#             rmhd_energy_solve.o

