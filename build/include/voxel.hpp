#ifndef __VOXEL__
#define __VOXEL__

struct voxel {
  double density;
  bool isAtSurface = false;
  bool isInSurface = false;
  // x, y, and z coordinates
  double x, y, z;
  // indices of the x, y, and z positions in the cube.  These are useful when
  // referencing other voxels by index (as in finding surface voxels)
  double xi, yi, zi;
};

bool compareVoxel(const voxel &a, const voxel &b);

#endif
