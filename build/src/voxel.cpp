#include "voxel.hpp"

/* Compare two voxels. */
bool compareVoxel(const voxel &a, const voxel &b) {
  return a.density < b.density;
}

/* Sum the densities in a vector of voxels */
double sumVoxelDensity(const std::vector<voxel> * voxels) {
  // accumulator
  double summ = 0.0;

  // sum the contents of the vector
  for (int ii = voxels.begin(); ii < voxels.end(); ii++) {
    summ += voxels[ii].density;
  }

  return summ;
}
