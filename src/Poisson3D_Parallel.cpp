#include "Poisson3D_Parallel.hpp"


void
Poisson3DParallel::write_csv(const long int elapsed_time, int iterations) const{
    std::ofstream csv_file;
    csv_file.open("../results/results.csv", std::ios_base::app);

    if (!csv_file.is_open()) {
        std::cerr<<"Error! Cannot open csv file! ";
        return ;
    }

    csv_file << extractFileName(mesh_file_name) << ","
             << p_value <<","
             << r <<","
             << preconditioner_name << ","
             << elapsed_time<<","
             << iterations <<","
             << symmetric <<","
             << mpi_size
             << std::endl;
    csv_file.close();

    pcout<< "CSV write successful."<<std::endl;
}

std::string
Poisson3DParallel::extractFileName(const std::string& filePath) {

    size_t lastSlash = filePath.find_last_of("/\\");


    std::string fileName = filePath.substr(lastSlash + 1);


    size_t lastDot = fileName.find_last_of(".");


    return fileName.substr(0, lastDot);
}

void
Poisson3DParallel::manage_flags(int argc, char **argv) {


    /*
     * DEFAULT INITIALIZATIONS OF PARAMETERS
     */

    std::string mesh_name_no_path = "mesh-cube-20.msh";
    preconditioner_name = "identity";
    p_value = 2;
    r = 1;
    symmetric = true;
    std::string user_choice_for_coefficient_symmetry;


    /*
     * MANAGING COMMAND LINE FLAGS
     */

    int opt;
    const char *options = "hp:m:r:P:s:";

    while ((opt = getopt(argc, argv, options))!=-1){
        switch (opt) {
            case 'h':
                pcout << "Help message "<<std::endl;
                break;

            case 'p':
                p_value = std::stod(optarg);
                pcout<<"Setting p="<<p_value<<std::endl;
                break;

            case 'm':
                mesh_name_no_path = optarg;
                pcout<<"Setting mesh="<<mesh_name_no_path<<std::endl;
                break;

            case 'P':
                preconditioner_name = optarg;
                pcout<<"Setting preconditioner=" <<preconditioner_name<<std::endl;
                break;

            case 'r':
                r = std::stoi(optarg);
                pcout<<"Setting r="<<r<<std::endl;
                break;

            case 's':
                user_choice_for_coefficient_symmetry = optarg;
                if (user_choice_for_coefficient_symmetry == "no") {
                    symmetric = false;
                    pcout<<"Initializing randomly an unsymmetric diffusion coefficient " << std::endl;
                }
                else
                    pcout<<"Initializing non-randomly a symmetric diffusion coefficient "<< std::endl;

                break;
            case '?':
                std::cerr << "Unknown command line option\n";
                return;

            default:
                pcout<<"default"<<std::endl;
                break;
        }

    }

    mesh_file_name = "../mesh/" + mesh_name_no_path;

    if(symmetric)
        initialize_diffusion_coefficient_symmetric(p_value);
    else
        initialize_diffusion_coefficient(p_value);

}


void
Poisson3DParallel::setup()
{
  pcout << "===============================================" << std::endl;

  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    // First we read the mesh from file into a serial (i.e. not parallel)
    // triangulation.
    Triangulation<dim> mesh_serial;

    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(mesh_serial);

      std::ifstream grid_in_file(mesh_file_name);
      grid_in.read_msh(grid_in_file);
    }

    // Then, we copy the triangulation into the parallel one.
    {
      GridTools::partition_triangulation(mpi_size, mesh_serial);
      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
      mesh.create_triangulation(construction_data);
    }

    // Notice that we write here the number of *global* active cells (across all
    // processes).
    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space. This is the same as in serial codes.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
   }
    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(r + 1);

    pcout << "  Quadrature points per boundary cell = "
              << quadrature_boundary->size() << std::endl;

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We retrieve the set of locally owned DoFs, which will be useful when
    // initializing linear algebra classes.
    locally_owned_dofs = dof_handler.locally_owned_dofs();

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // To initialize the sparsity pattern, we use Trilinos' class, that manages
    // some of the inter-process communication.
    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);

    // After initialization, we need to call compress, so that all process
    // retrieve the information they need for the rows they own (i.e. the rows
    // corresponding to locally owned DoFs).
    sparsity.compress();

    // Then, we use the sparsity pattern to initialize the system matrix. Since
    // the sparsity pattern is partitioned by row, so will the matrix.
    pcout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity);

    // Finally, we initialize the right-hand side and solution vectors.
    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  }
}

void
Poisson3DParallel::assemble()
{
  pcout << "===============================================" << std::endl;

  pcout << "  Assembling the linear system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);


  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // If current cell is not owned locally, we skip it.
      if (!cell->is_locally_owned())
        continue;

      // On all other cells (which are owned by current process), we perform the
      // assembly as usual.

      fe_values.reinit(cell);

      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Diffusion term.
                  cell_matrix(i, j) += diffusion_coefficient.value(
                                         fe_values.quadrature_point(q)) // mu(x)
                                       * fe_values.shape_grad(i, q)     // (I)
                                       * fe_values.shape_grad(j, q)     // (II)
                                       * fe_values.JxW(q);              // (III)
                }

              cell_rhs(i) += forcing_term.value(fe_values.quadrature_point(q)) *
                             fe_values.shape_value(i, q) * fe_values.JxW(q);
            }
        }
        // If the cell is adjacent to the boundary...
      

      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // Each process might have written to some rows it does not own (for instance,
  // if it owns elements that are adjacent to elements owned by some other
  // process). Therefore, at the end of the assembly, processes need to exchange
  // information: the compress method allows to do this.
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  // Boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;

    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    for (unsigned int i = 0; i < 6; ++i)
      boundary_functions[i] = &function_g;

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, true);
  }
}

void
Poisson3DParallel::solve()
{
  pcout << "===============================================" << std::endl;

  SolverControl solver_control(10000, 1e-6 * system_rhs.l2_norm());

  // The linear solver is basically the same as in serial, in terms of
  // interface: we only have to use appropriate classes, compatible with
  // Trilinos linear algebra.
  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
  SolverGMRES<TrilinosWrappers::MPI::Vector> GMRESsolver(solver_control);

  pcout << "  Solving the linear system" << std::endl;
  const auto t0 = std::chrono::high_resolution_clock::now();




    if (preconditioner_name == "identity")
    {
        std::cout<<"Using preconditioner identity"<<std::endl;
        TrilinosWrappers::PreconditionIdentity preconditioner;
        solver.solve(system_matrix, solution, system_rhs, preconditioner);
    }else if (preconditioner_name == "jacobi"){
        std::cout<<"Using preconditioner jacobi"<<std::endl;
        TrilinosWrappers::PreconditionJacobi preconditioner;
        preconditioner.initialize(system_matrix);
        solver.solve(system_matrix, solution, system_rhs, preconditioner);
    }else if (preconditioner_name == "ssor"){
        std::cout<<"Using preconditioner ssor"<<std::endl;
        TrilinosWrappers::PreconditionSSOR preconditioner;
        preconditioner.initialize(
                system_matrix,
                TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0)
        );
        solver.solve(system_matrix, solution, system_rhs, preconditioner);

    }else if (preconditioner_name == "sor"){

        std::cout<<"Using preconditioner sor"<<std::endl;
        std::cout<<"Not symmetric preconditioner, hence solving with GMRES and not GC"<<std::endl;
        TrilinosWrappers::PreconditionSOR preconditioner;
        preconditioner.initialize(
                system_matrix, TrilinosWrappers::PreconditionSOR::AdditionalData(1.0)
        );
        GMRESsolver.solve(system_matrix, solution, system_rhs, preconditioner);
    }
    else{
        std::cerr<<"Error! Preconditioner not supported!"<<std::endl;
        std::exit(-1);
    }

    const auto t1= std::chrono::high_resolution_clock::now();


    const auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    pcout<<"Elapsed time for solve phase: "
             << dt << " [ms] " << std::endl;

    pcout << "  " << solver_control.last_step() << "  iterations" << std::endl;

  if (mpi_rank == 0)
      write_csv(dt, solver_control.last_step());

}

void
Poisson3DParallel::output() const
{
  pcout << "===============================================" << std::endl;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // To correctly export the solution, each process needs to know the solution
  // DoFs it owns, and the ones corresponding to elements adjacent to the ones
  // it owns (the locally relevant DoFs, or ghosts). We create a vector to store
  // them.
  TrilinosWrappers::MPI::Vector solution_ghost(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               MPI_COMM_WORLD);

  // This performs the necessary communication so that the locally relevant DoFs
  // are received from other processes and stored inside solution_ghost.
  solution_ghost = solution;

  // Then, we build and fill the DataOut class as usual.
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution_ghost, "solution");

  // We also add a vector to represent the parallel partitioning of the mesh.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string output_file_name = "output-" + mesh_path.stem().string();

  // Finally, we need to write in a format that supports parallel output. This
  // can be achieved in multiple ways (e.g. XDMF/H5). We choose VTU/PVTU files,
  // because the interface is nice and it is quite robust.
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      0,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;

  pcout << "===============================================" << std::endl;
}