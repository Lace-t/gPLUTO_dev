#ifndef GPLUTO_GRID_SPEC_HPP
#define GPLUTO_GRID_SPEC_HPP
// POD layout: Runtime is broadcast as bytes in MPI builds.
constexpr int GRID_MAX_PATCHES = 16;
struct GridSpec {
  int count = 0; // zero selects legacy three-field syntax / MappedGrid()
  int cells[GRID_MAX_PATCHES] = {};
  char type[GRID_MAX_PATCHES] = {}; // u, s, +, -
  double edge[GRID_MAX_PATCHES+1] = {};
};
#endif
