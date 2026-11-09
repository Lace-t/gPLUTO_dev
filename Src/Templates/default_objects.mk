# Default object and header files for gPLUTO

DEFAULT_OBJ = \
  arrays.o array_reconstruct.o bin_io.o boundary.o cmd_line_opt.o \
  data2sweep.o debug_tools.o euler_step.o \
  fallback.o flag_boundary.o flag_exchange.o flag_shock.o flat_states.o \
  get_nghost.o globals.o \
  init.o initialize.o input_data.o internal_boundary.o jet_domain.o \
  main.o mapped_grid_tools.o \
  mappers3D.o mean_mole_weight.o mem_usage.o next_time_step.o \
  nvtx.o output_grid.o output_log.o parallel_exchange_async.o \
  parallel_exchange_sync.o parse_file.o rbox.o   \
  ring_average.o ring_average_rec.o              \
  restart.o rotate.o runtime_setup.o source.o    \
  set_geometry.o set_grid.o set_output.o set_var_names.o \
  show_config.o split_source.o startup.o sweep_alloc.o \
  timing.o tools.o write_data.o write_tab.o write_vtk.o \
  vector_curl.o userdef_output.o

DEFAULT_HEADERS = Array.hpp definitions.hpp eos.hpp gpluto.hpp \
                  macros.hpp metric.hpp plm_coeffs.hpp prototypes.hpp \
                  states.hpp structs.hpp

