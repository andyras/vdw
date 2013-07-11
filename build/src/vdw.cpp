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
#include "params.hpp"

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
  PARAMETERS params;	// struct of params

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
  params.nx = nx;
  params.ny = ny;
  params.nz = nz;

  // If the sign of the number of voxels is positive, then the units are
  // Bohr, if negative then Angstrom.
  // TODO: Gaussian manual contradicts this code; in the G09 manual (nx < 0)
  // means Bohr.
  // INTERNALLY THIS PROGRAM USES ANGSTROM
  if (nx < 0) {
    ANGSTROM = true;
    nx = abs(nx);
    ny = abs(ny);
    nz = abs(nz);
  }
  // convert to Angstrom
  else {
    for (int ii = 0; ii < 3; ii++) {
      origin[ii] *= a0;
      xv[ii] *= a0;
      yv[ii] *= a0;
      zv[ii] *= a0;
    }
  }
  std::cout << "The origin is at     " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
  std::cout << "The x unit vector is " << xv[0] << " " << xv[1] << " " << xv[2] << std::endl;
  std::cout << "The y unit vector is " << yv[0] << " " << yv[1] << " " << yv[2] << std::endl;
  std::cout << "The z unit vector is " << zv[0] << " " << zv[1] << " " << zv[2] << std::endl;

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
  params.lx = lx;
  params.ly = ly;
  params.lz = lz;
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
    // shift all atomic coordinates
    // TODO why are these being shifted the wrong way?
    atoms[ii].x -= origin[0];
    atoms[ii].y -= origin[1];
    atoms[ii].z -= origin[2];
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
	voxels[idx].density = fabs(scratch);	// electron density
	voxels[idx].xi = ii;			// x, y, z indices
	voxels[idx].yi = jj;
	voxels[idx].zi = kk;
	voxels[idx].x = (origin[0] + (ii-1)*xv[0] + (jj-1)*xv[1] + (kk-1)*xv[2]);	// x, y, z coordinates
	voxels[idx].y = (origin[0] + (ii-1)*yv[0] + (jj-1)*yv[1] + (kk-1)*yv[2]);
	voxels[idx].z = (origin[0] + (ii-1)*zv[0] + (jj-1)*zv[1] + (kk-1)*zv[2]);
      }
    }
  }
  cube.close();

  std::cout << "Done reading " << cubeFileName << "." << std::endl;

  //// sum vector of voxels
  std::cout << "Summing total electron density..." << std::endl;
  double totalDensity = sumVoxelDensity(voxels);
  std::cout << "Total electron density is " << totalDensity << std::endl;
  std::cout << "Total electron population is " << totalDensity*vol << std::endl;

  //// calculate cutoff density (fraction of total density)
  double cutoffDensity = densityCutoffPercent*totalDensity;
  std::cout << "Cutoff electron density is " << cutoffDensity << std::endl;

  //// sort vector of voxels
  // first copy voxel vector, for reference to find indices
  std::vector<voxel> voxelsCopy (voxels);

  // highest density values will now be at the beginning of voxels vector
  std::sort(voxels.begin(), voxels.end(), compareVoxel);

  //// find voxels which are within the isosurface
  // starting with highest density, accumulate until cutoff is reached
  double sumDensity = 0.0;
  // surfaceIndex is a variable which will define the range of voxels in the
  // surface surrounding the cutoff density.
  // At the end of the loop, surfaceIndex will have a value that can be used in
  // a loop, e.g. 'for (int ii = 0; ii < surfaceIndex; ii++)' will loop over
  // elements within the overall surface.
  int surfaceIndex = 0;
  while (sumDensity < cutoffDensity) {
    sumDensity += voxels[surfaceIndex].density;
    voxels[surfaceIndex].isInSurface = true;
    //std::cout << "ii " << surfaceIndex << " density " << voxels[surfaceIndex].density << " sumDensity " << sumDensity << std::endl;
    surfaceIndex++;
  }
  //for (surfaceIndex = 0; sumDensity < cutoffDensity; sumDensity += voxels[surfaceIndex++].density) {
    //voxels[surfaceIndex].isInSurface = true;
  //}
  // find the closest possible value of the cutoff density
  if (fabs(sumDensity - cutoffDensity) > fabs(sumDensity - voxels[surfaceIndex-1].density - cutoffDensity)) {
    sumDensity -= voxels[surfaceIndex-1].density;
    surfaceIndex -= 1;
  }
  std::cout << "Sum of voxels approximating cutoff density: " << sumDensity << std::endl;
  std::cout << "Number of voxels within cutoff density: " << surfaceIndex << std::endl;
  std::cout << "Volume within cutoff density: " << vol*surfaceIndex << std::endl;

  //// find electron density on the atom of interest
  std::cout << "Finding voxels closest to atom of interest" << std::endl;
  int voxelsInAtom = surfaceIndex;	// count of the voxels which are in the atom of interest
  double voxelAOIDistance;		// distance from a voxel to the atom of interest
  for (int ii = 0; ii < surfaceIndex; ii++) {
    for (int jj = 0; jj < natoms; jj++) {
      voxelAOIDistance = voxelAtomDistance(&voxels[ii], &atoms[aoi]);
      //if (ii <= 10) {
	//std::cout << "Distance from voxel " << ii << " to atom " << aoi << ": " << voxelAOIDistance << std::endl;
	//std::cout << "Offset is "
	/*
      std::cout << (voxels[ii].x - atoms[aoi].x)
		  << " " << (voxels[ii].y - atoms[aoi].y)
		  << " " << (voxels[ii].z - atoms[aoi].z)
		  << " " << voxels[ii].density << std::endl;
		  */
      //}
      if ((jj != aoi) && (voxelAtomDistance(&voxels[ii], &atoms[jj]) < voxelAOIDistance)) {
	voxels[ii].isInAtom = false;
	voxelsInAtom -= 1;
	break;
      }
    }
  }
  std::cout << "Number of voxels within atom of interest: " << voxelsInAtom << std::endl;

  //// find voxels at the surface of the atom of interest
  std::vector<voxel> surfaceVoxels;
  for (int ii = 0; ii < surfaceIndex; ii++) {
    if (checkIfSurfaceVoxel(&voxels[ii], voxelsCopy, &params)) {
      voxels[ii].isAtSurface = true;
      surfaceVoxels.push_back(voxels[ii]);
    }
  }
  //// find distances between voxels on surface of atom
  //// find largest distance between voxels on atoms
  //// write new .cube file, just with density on atom
  return 0;
}
