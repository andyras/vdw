#include "voxel.hpp"

/* Compare two voxels. */
bool compareVoxel(const voxel &a, const voxel &b) {
  return a.density < b.density;
}
