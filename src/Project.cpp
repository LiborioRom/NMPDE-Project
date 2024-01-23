#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson3D.hpp"

// Main function.
int
main(int argc, char * argv[])
{

  Poisson3D problem;
  problem.manage_flags(argc, argv);


  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}