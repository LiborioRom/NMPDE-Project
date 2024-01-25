#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include <mpi.h>

#include "Poisson3D.hpp"

// Main function.
int
main(int argc, char * argv[])
{
  MPI_Init(&argc, &argv);

   int world_size, world_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Only let rank 0 execute the MPI code
  if (world_rank == 0) {
      Poisson3D problem;
      problem.manage_flags(argc, argv);
      problem.setup();
      problem.assemble();
      problem.solve();
      problem.output();
    //  std::cout << "MPI code executed by rank " << world_rank << std::endl;
  }


 MPI_Finalize();
  return 0;
}

/*
mpiexec -n 4 ./Project_parallel -p 4 -P amg -r 2 -m mesh-cube-40.msh -s 0
./Project -p 4 -P amg -r 2 -m mesh-cube-40.msh -s 0

-p : set p value (10^p)
-P : set preconditionner
-r : set degree
-m : set mesh file (without path)
-s : Choose between symmetric (any value) /unsymmetric ("no")

relation between continuity difussion problem coefficiants and discrete problem
*/













/*
we need to add these...
-ilu 
-ilut 
-amg
-blockwise_direct


poisson3d serial :

PreconditionIdentity
PreconditionJacobi
PreconditionSOR
PreconditionSSOR

PreconditionAMG (to be checked...)



poisson3d parallel :

PreconditionIdentity
PreconditionJacobi
PreconditionSOR
PreconditionSSOR
*/