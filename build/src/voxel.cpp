#include "voxel.hpp"

//#define DEBUG_VOXEL

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

#ifdef DEBUG_VOXEL
  std::cout << "Voxel is " << std::setw(4) << v->xi << " "
                           << std::setw(4) << v->yi << " "
                           << std::setw(4) << v->zi << std::endl;
#endif
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
#ifdef DEBUG_VOXEL
	  std::cout << "Neighbor " << std::setw(4) << ii << " "
				   << std::setw(4) << jj << " "
				   << std::setw(4) << kk << " is ";
#endif
	  if (!(vs[ii*p->ny*p->nz + jj*p->nz + kk].isInAtom)) {
#ifdef DEBUG_VOXEL
	    std::cout << "NOT in atom" << std::endl;
#endif
	    return true;
	  }
#ifdef DEBUG_VOXEL
	  std::cout << "in atom" << std::endl;
#endif
	}
      }
    }
  }
  // All other cases mean the voxel is not at the surface.
  return false;
}

/* Finds the distance between two voxels. */
double voxelDistance(const struct voxel * v1, const struct voxel * v2) {
  return sqrt(pow(v1->x - v2->x, 2) + pow(v1->y - v2->y, 2) + pow(v1->z - v2->z, 2));
}

/* Decides whether a voxel is in an atom of interest. */
bool voxelInAtom(const struct voxel * v, const int a, std::vector<atom> &as) {
    bool isInAtom = true;
    double voxelAOIDistance = 1e20;
    for (int jj = 0; jj < as.size(); jj++) {
      voxelAOIDistance = voxelAtomDistance(v, &as[a]);
      // a voxel is in the atom if no other atoms are closer
      if ((jj != a) && (voxelAtomDistance(v, &as[jj]) < voxelAOIDistance)) {
	isInAtom = false;
	break;
      }
    }
  return isInAtom;
}
