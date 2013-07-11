#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "voxel.hpp"
#include "atom.hpp"
#include "constants.hpp"

int main(int argc, char ** argv) {
  bool printHelp = false;
  // variable for parsing inputs
  int c;

  // print help if no argument given
  if (argc == 1) {
    printHelp = true;
  }

  int aoi;			// index of the atom of interest
  double densityCutoffPercent;	// density cutoff value (fraction of total electron density)

  //// parse input parameters
  while ((c = getopt(argc, argv, "a:d:h")) != -1) {
    switch (c) {
      case 'a':
	aoi = atoi(optarg);
	std::cout << "Atom index is " << aoi << std::endl;
	break;
      case 'd':
	densityCutoffPercent = atof(optarg);
	std::cout << "Density cutoff is " << densityCutoffPercent << std::endl;
	break;
      case 'h':
	printHelp = true;
	break;
      case '?':
	if (optopt == 'a') {
	  fprintf(stderr, "Option -%c requires an atom index.\n", optopt);
	}
	if (optopt == 'd') {
	  fprintf(stderr, "Option -%c requires a density cutoff (0.0 < cutoff < 1.0).\n", optopt);
	}
	else if (isprint(optopt)) {
	  fprintf(stderr, "Unknown option -%c.\n", optopt);
	}
	else {
	  fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
	}
	return 1;
	default:
	  continue;
    }
  }

  if (printHelp) {
    std::cout << std::endl;
    std::cout << "vdw: A program for calculating some van der Waals properties of a molecule." << std::endl;
    std::cout << "     It finds the van der Waals radius of an atom in the molecule." << std::endl;
    std::cout << "     Much of the code (especially .cube parsing) was borrowed from" << std::endl;
    std::cout << "     Nathaniel Swenson's 'bq' code for estimating vdW volumes." << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Usage: vdw [-h] -a <atom index> -d <density cutoff> <cube file>" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Only the help flag is optional. A cube file is required as input." << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Flags:" << std::endl;
    std::cout << "-a: the index (starting at 0) of the atom you wish to examine" << std::endl;
    std::cout << "-d: the density cutoff (between 0 and 1)" << std::endl;
    std::cout << "-h: print this help message" << std::endl;
    std::cout << "" << std::endl;
    return 2;
  }

  if (argc == optind) {
    std::cerr << "ERROR: cube file must be specified.\n";
    return 1;
  }

  // get the input .cube file name
  char * cubeFileName = argv[optind];
  std::cout << "The input filename is '" << cubeFileName << "'" << std::endl;

  //// read in .cube file
  std::ifstream cube;
  std::string cl;	// current line
  int natoms;		// number of atoms
  int nx, ny, nz;	// number of voxels in each dimension
  double origin [3];	// coordinates of origin
  double xv [3];	// vectors of x, y, z directions
  double yv [3];
  double zv [3];
  double lx, ly, lz;	// length of voxel vector in each dimension
  double vol;		// voxel volume
  double vol_a, vol_b;	// not sure what these are
  bool ANGSTROM = false;// switch for Angstrom or Bohr units of length
  double UNITS = a0;	// units of length

  std::cout << "Reading input file " << cubeFileName << "..." << std::endl;
  // open the file
  cube.open(cubeFileName, std::ios::in);

  // the first two lines are comments; ignore them
  std::getline(cube,cl);
  std::getline(cube,cl);

  // Line 3 is the number of atoms followed by the origin of the volumetric data.
  cube >> natoms >> origin[0] >> origin[1] >> origin[2];
  if (aoi > (natoms + 1)) {
    std::cerr << "ERROR: index of atom of interest is greater than the number of atoms." << std::endl;
    return 1;
  }

  // Lines 4-6 give the number of voxels along each axis (x,y,z) followed
  // by the axis vector.
  cube >> nx >> xv[0] >> xv[1] >> xv[2];
  cube >> ny >> yv[0] >> yv[1] >> yv[2];
  cube >> nz >> zv[0] >> zv[1] >> zv[2];
  std::cout << "Grid dimensions are " << nx << " x " << ny << " x " << nz << " points." << std::endl;

  // If the sign of the number of voxels is positive, then the units are
  // Bohr, if negative then Angstrom.
  // TODO: Gaussian manual contradicts this code; in the G09 manual (nx < 0)
  // means Bohr.
  // INTERNALLY THIS PROGRAM USES ANGSTROM
  if (nx < 0) {
    ANGSTROM = true;
    UNITS = 1.0;
    nx = abs(nx);
    ny = abs(ny);
    nz = abs(nz);
  }

  // Compute the voxel volume
  lx = sqrt(pow(xv[0], 2) + pow(xv[1], 2) + pow(xv[2], 2));
  ly = sqrt(pow(yv[0], 2) + pow(yv[1], 2) + pow(yv[2], 2));
  lz = sqrt(pow(zv[0], 2) + pow(zv[1], 2) + pow(zv[2], 2));
  // if the lengths are in Bohr, convert to Angstrom
  if (!ANGSTROM) {
    lx *= a0;
    ly *= a0;
    lz *= a0;
  }
  vol = lx*ly*lz;
  std::cout << "Grid extends over " << nx*lx << " x " << ny*ly << " x " << nz*lz << " Angstrom." << std::endl;
  std::cout << "Grid extends over " << nx*lx/a0 << " x " << ny*ly/a0 << " x " << nz*lz/a0 << " Bohr." << std::endl;
  std::cout << "Voxel volume is " << lx*ly*lz << " Angstrom^3 (" << lx*ly*lz/pow(a0,3) << " Bohr^3)." << std::endl;

  // TODO what does this do
  if (ANGSTROM) {
    vol_a = vol;
    vol_b = vol/(pow(a0,3));
  }
  else {
    vol_a = vol*pow(a0,3);
    vol_b = vol_a;
  }

  // The remaining lines in the header are the atomic coordinates.
  std::vector<atom> atoms;
  atoms.resize(natoms);
  double scratch;		// variable to capture junk
  for (int ii = 0; ii < natoms; ii++) {
    cube >> atoms[ii].atomicNo >> scratch >> atoms[ii].x >> atoms[ii].y >> atoms[ii].z;
    if (!ANGSTROM) {
      atoms[ii].x *= a0;
      atoms[ii].y *= a0;
      atoms[ii].z *= a0;
    }
  }
  std::cout << std::endl;
  std::cout << "Coordinates of atoms:" << std::endl;
  for (int ii = 0; ii < atoms.size(); ii++) {
    std::cout << atoms[ii].atomicNo;
    std::cout << " " << std::setw(8) << std::scientific << atoms[ii].x;
    std::cout << " " << std::setw(8) << std::scientific << atoms[ii].y;
    std::cout << " " << std::setw(8) << std::scientific << atoms[ii].z;
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // The rest of the file is the volumetric data.
  std::vector<voxel> voxels;
  voxels.resize(nx*ny*nz);
  int idx;			// index in voxels array

  for (int ii = 0; ii < nx; ii++) {
    for (int jj = 0; jj < ny; jj++) {
      for (int kk = 0; kk < nz; kk++) {
	idx = ii*ny*nz + jj*nz + kk;
	cube >> scratch;
	voxels[idx].density = vol*scratch;		// electron density times volume element
	voxels[idx].xi = ii;				// x, y, z indices
	voxels[idx].yi = jj;
	voxels[idx].zi = kk;
	voxels[idx].x = (origin[0] + ii*lx)*UNITS;	// x, y, z coordinates
	voxels[idx].y = (origin[1] + jj*ly)*UNITS;
	voxels[idx].z = (origin[2] + kk*lz)*UNITS;
      }
    }
  }
  cube.close();

  std::cout << "Done reading " << cubeFileName << "." << std::endl;

  //// sum vector of voxels
  std::cout << "Summing total electron density..." << std::endl;
  double totalDensity = sumVoxelDensity(voxels);
  std::cout << "Total electron density is " << totalDensity << std::endl;

  //// calculate cutoff density (fraction of total density)
  double cutoffDensity = densityCutoffPercent*totalDensity;

  //// sort vector of voxels
  // first copy voxel vector, for reference to find indices
  std::vector<voxel> voxelsCopy (voxels);

  // highest density values will now be at the beginning of voxels vector
  std::sort(voxels.begin(), voxels.end(), compareVoxel);

  //// find voxels which are within the isosurface
  // starting with highest density, accumulate until cutoff is reached

  //// find electron density on the atom of interest
  //// find voxels at the surface of the atom of interest
  //// find distances between voxels on surface of atom
  //// find largest distance between voxels on atoms
  return 0;
}
