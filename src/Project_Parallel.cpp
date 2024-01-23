#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson3D_Parallel.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  // This object calls MPI_Init when it is constructed, and MPI_Finalize when it
  // is destroyed. It also initializes several other libraries bundled with
  // dealii (e.g. p4est, PETSc, ...).

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);


  Poisson3DParallel problem;
  problem.manage_flags(argc, argv);
  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}