#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "voxel.hpp"
#include "atom.hpp"

int main(int argc, char ** argv) {
  bool printHelp = false;
  // variable for parsing inputs
  int c;

  // print help if no argument given
  if (argc == 1) {
    printHelp = true;
  }

  //// parse input parameters
  while ((c = getopt(argc, argv, "a:d:h")) != -1) {
    switch (c) {
      case 'a':
	std::cout << "Atom index is " << optarg << std::endl;
	break;
      case 'd':
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

  //// read in .cube file
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
