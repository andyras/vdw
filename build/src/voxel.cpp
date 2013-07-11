#include "voxel.hpp"

/* Compare two voxels.
 * When used with std::sort, gives list sorted with biggest first.
 */
bool compareVoxel(const voxel &a, const voxel &b) {
  return a.density > b.density;
}

/* Sum the densities in a vector of voxels */
double sumVoxelDensity(std::vector<voxel> &voxels) {
  // accumulator
  double summ = 0.0;

  // sum the contents of the vector
  for (int ii = 0; ii < voxels.size(); ii++) {
    summ += voxels[ii].density;
  }

  return summ;
}

/* Finds the distance between a voxel and an atom. */
double voxelAtomDistance(const struct voxel * v, const struct atom * a) {
  return sqrt(pow(v->x - a->x, 2) + pow(v->y - a->y, 2) + pow(v->z - a->z, 2));
}

/* Decides whether voxel is at the surface or not. */
bool checkIfSurfaceVoxel(const struct voxel * v, std::vector<voxel> &vs,
                         const struct PARAMETERS * p) {
  // Check the surrounding voxels to see if they are also in the atom.
  // If any are not, then the voxel is at the surface.
  // If any do not exist (i.e. the voxel is at the edge of the cube),
  // then the voxel is at the surface.

  std::cout << "Voxel is " << std::setw(4) << v->xi << " "
                           << std::setw(4) << v->yi << " "
                           << std::setw(4) << v->zi << std::endl;
  // Check first whether voxel is at edge of box
  if (   (v->xi == 0) || (v->xi == (p->nx-1))
      || (v->yi == 0) || (v->yi == (p->ny-1))
      || (v->zi == 0) || (v->zi == (p->nz-1))) {
    return true;
  }
  // Next check if neighbors are not in the atom
  else {
    for (int ii = (v->xi-1); ii < (v->xi+2); ii++) {
      for (int jj = (v->yi-1); jj < (v->yi+2); jj++) {
	for (int kk = (v->zi-1); kk < (v->zi+2); kk++) {
	  std::cout << "Neighbor " << std::setw(4) << ii << " "
				   << std::setw(4) << jj << " "
				   << std::setw(4) << kk << " is ";
	  if (!(vs[ii*p->ny*p->nz + jj*p->nz + kk].isInAtom)) {
	    std::cout << "NOT in atom" << std::endl;
	    return true;
	  }
	  std::cout << "in atom" << std::endl;
	}
      }
    }
  }
  // All other cases mean the voxel is not at the surface.
  return false;
}
