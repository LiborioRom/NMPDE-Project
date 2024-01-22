#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson3D.hpp"

// Main function.
int
main(int argc, char * argv[])
{

  std::string mesh_name = argv[1];
  const std::string mesh_file_name =
    "../mesh/" + mesh_name;   //paralepiped.geo.msh


  const unsigned int r = 1;

  Poisson3D problem(mesh_file_name, r);
  problem.manage_flags(argc, argv);


  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}