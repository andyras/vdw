#ifndef __VOXEL__
#define __VOXEL__

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "atom.hpp"
#include "params.hpp"

struct voxel {
  double density;
  // voxel is in isosurface
  bool isInIsoSurface = false;
  // voxel is in atom of interest
  bool isInAtom = true;
  // voxel is at surface of atom of interest
  bool isAtSurface = false;
  // x, y, and z coordinates
  double x, y, z;
  // indices of the x, y, and z positions in the cube.  These are useful when
  // referencing other voxels by index (as in finding surface voxels)
  int xi, yi, zi;
};

bool compareVoxel(const voxel &a, const voxel &b);

double sumVoxelDensity(std::vector<voxel> &voxels);

double voxelAtomDistance(const struct voxel * v, const struct atom * a);

bool checkIfSurfaceVoxel(const struct voxel * v, std::vector<voxel> &vs,
                         const struct PARAMETERS * p);

#endif
