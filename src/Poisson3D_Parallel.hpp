#ifndef POISSON_3D_HPP
#define POISSON_3D_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <Epetra_MpiComm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Export.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_LinearProblem.h>


#include <filesystem>
#include <fstream>
#include <iostream>

#include <vector>
#include <cmath> // For pow and sqrt functions
#include <random>


using namespace dealii;

/**
 * Class managing the differential problem.
 */
class Poisson3DParallel
{
public:

  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Diffusion coefficient.
  template <int dim>
  class DiffusionCoefficient : public Function<dim>
  {
  public:
      // Define a structure to hold sphere data
      struct Sphere
      {
          Point<dim> center;
          double radius;
      };
  
      // Constructor
      DiffusionCoefficient(const std::vector<Sphere>& spheres, double p)
      : spheres(spheres), p_value(p) {}

      DiffusionCoefficient(){}

      // Evaluation
      virtual double value(const Point<dim>& p, const unsigned int /*component*/ = 0) const override
      {
          for (const auto& sphere : spheres)
          {
              if (is_point_inside_sphere(p, sphere))
                  return std::pow(10.0, p_value);
          }
          return 1.0;
      }

  private:
      std::vector<Sphere> spheres;
      double p_value;

      bool is_point_inside_sphere(const Point<dim>& point, const Sphere& sphere) const
      {
          double distance_squared = 0.0;
          for (unsigned int i = 0; i < dim; ++i)
          {
              distance_squared += std::pow(point[i] - sphere.center[i], 2);
          }
          return distance_squared <= std::pow(sphere.radius, 2);
      }
  };

  // Forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    // Constructor.
    ForcingTerm()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };

    // Neumann boundary conditions.
    class FunctionH : public Function<dim>
    {
    public:
        // Constructor.
        FunctionH()
        {}

        // Evaluation:
        virtual double
        value(const Point<dim> /*&p*/, const unsigned int /*component*/ = 0) const
        {
            return 0.0;
        }
    };
// Dirichlet boundary conditions.
    class FunctionG : public Function<dim>
    {
    public:
        // Constructor.
        FunctionG()
        {}

        // Evaluation.
        virtual double
        value(const Point<dim> & /*p*/,
              const unsigned int /*component*/ = 0) const override
        {
            return 0.0;
        }
    };


    // Constructor.
  Poisson3DParallel():
    mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
    mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
    mesh(MPI_COMM_WORLD),
    pcout(std::cout, mpi_rank == 0)
  {};


  void initialize_diffusion_coefficient(double p_value) {
        // Crear un vector de esferas
        std::vector<DiffusionCoefficient<dim>::Sphere> spheres;

        // Generar esferas aleatorias
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis_center(0.0, 1.0);
        std::uniform_real_distribution<> dis_radius(0.01, 0.1); // Radios pequeños para asegurarse de que caben en el cubo

        for (int i = 0; i < 10; ++i) {
            Point<dim> center = {dis_center(gen), dis_center(gen), dis_center(gen)};
            double radius = dis_radius(gen);
            spheres.push_back({center, radius});
        }


        // Inicializar el coeficiente de difusión
        diffusion_coefficient = DiffusionCoefficient<dim>(spheres, p_value);
    };

  void initialize_diffusion_coefficient_symmetric(double p_value) {
        // Create a vector of spheres
        std::vector<DiffusionCoefficient<dim>::Sphere> spheres;
        int n=2;
        // Define symmetrically distributed centers
        double step = 1.0 / (n+1);  // 'n' is the number of steps to divide each dimension
        double radius = 0.05;       // Fixed radius for each sphere

        // Generate symmetrically placed spheres
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                for (int k = 1; k <= n; ++k) {
                    if (spheres.size() < 10) {
                        Point<dim> center = {i * step, j * step, k * step};
                        spheres.push_back({center, radius});
                    }
                }
            }
        }
      // Initialize the diffusion coefficient
      diffusion_coefficient = DiffusionCoefficient<dim>(spheres, p_value);
  };
    
  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve();

  // Output.
  void
  output() const;

  // Export to matrix market.
  void exportToMatrixMarket(const TrilinosWrappers::SparseMatrix &matrix, const std::string &filename);

  /*
   * CUSTOM PUBLIC MEMBERS (SOME CAN BE MADE PROTECTED)
   */

  void
  manage_flags(int argc, char ** argv);

  void
  write_csv(const long int elapsed_time, int iterations, double condition_number) const;

  static std::string
  extractFileName(const std::string& filePath);

protected:
  // Path to the mesh file.
  std::string mesh_file_name;

  /*
   * CUSTOM PROTECTED MEMBERS
   */

  std::string preconditioner_name;

  bool symmetric;

  double p_value;

  // Polynomial degree.
  unsigned int r; //Modified from const unsigned

  /********END CUSTOM MEMBERS*******/

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Diffusion coefficient.
  DiffusionCoefficient<dim> diffusion_coefficient;

  // Forcing term.
  ForcingTerm forcing_term;

    // h(x).
    FunctionH function_h;

    // g(x).
    FunctionG function_g;


  // Triangulation. The parallel::fullydistributed::Triangulation class manages
  // a triangulation that is completely distributed (i.e. each process only
  // knows about the elements it owns and its ghost elements).
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // System matrix.
  TrilinosWrappers::SparseMatrix system_matrix;

  // System right-hand side.
  TrilinosWrappers::MPI::Vector system_rhs;

  // System solution.
  TrilinosWrappers::MPI::Vector solution;

  // Parallel output stream.
  ConditionalOStream pcout;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;
};

#endif