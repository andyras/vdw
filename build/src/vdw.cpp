#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <iostream>
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

  int aoi;		// index of the atom of interest
  double densityCutoff;	// density cutoff value (fraction of total electron density)

  //// parse input parameters
  while ((c = getopt(argc, argv, "a:d:h")) != -1) {
    switch (c) {
      case 'a':
	aoi = atoi(optarg);
	std::cout << "Atom index is " << aoi << std::endl;
	break;
      case 'd':
	densityCutoff = atof(optarg);
	std::cout << "Density cutoff is " << optarg << std::endl;
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

  // Compute the voxel volume
  lx = sqrt(pow(xv[0], 2) + pow(xv[1], 2) + pow(xv[2], 2));
  ly = sqrt(pow(yv[0], 2) + pow(yv[1], 2) + pow(yv[2], 2));
  lz = sqrt(pow(zv[0], 2) + pow(zv[1], 2) + pow(zv[2], 2));
  vol = lx*ly*lz;

  // If the sign of the number of voxels is positive, then the units are
  // Bohr, if negative then Angstrom.
  if (nx < 0) {
    ANGSTROM = true;
    UNITS = 1.0;
    nx *= -1;
    ny *= -1;
    nz *= -1;
  }

  // TODO what does this do
  if (ANGSTROM) {
    vol_a = vol;
    vol_b = vol/(pow(a0,3));
  }
  else {
    vol_a = vol*pow(a0,3);
    vol_b = vol_a;
  }

  //// sum vector of voxels
  //// calculate cutoff density (fraction of total density)
  //// sort vector of voxels
  //// find voxels which are within the isosurface
  //// find electron density on the atom of interest
  //// find voxels at the surface of the atom of interest
  //// find distances between voxels on surface of atom
  //// find largest distance between voxels on atoms
  return 0;
}
