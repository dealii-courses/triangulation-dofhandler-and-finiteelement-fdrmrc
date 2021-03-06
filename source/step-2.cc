/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-2
 */


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <fstream>

using namespace dealii;


void
show_statistics(const SparsityPattern &sparsity_pattern)
{
  Vector<int> v;
  double      sum = 0.0;
  for (unsigned int i = 0; i < sparsity_pattern.n_cols(); ++i)
    {
      sum += sparsity_pattern.row_length(i);
    }

  std::cout << "Statistics for the current sparsity pattern:"
            << "\n"
            << "bandwidth: " << sparsity_pattern.bandwidth() << "\t"
            << " number of unknowns : " << sparsity_pattern.n_cols() << "\t"
            << "average number of of entries per row"
            << sum / sparsity_pattern.n_cols() << "\n";
}

void
make_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 5);

  static const SphericalManifold<2> manifold_description(center);
  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold_description);

  for (unsigned int step = 0; step < 3; ++step)
    {
      Triangulation<2>::active_cell_iterator cell =
                                               triangulation.begin_active(),
                                             endc = triangulation.end();

      for (; cell != endc; ++cell)
        for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
          {
            const double distance_from_center =
              center.distance(cell->vertex(v));

            if (std::fabs(distance_from_center - inner_radius) < 1e-10)
              {
                cell->set_refine_flag();
                break;
              }
          }

      triangulation.execute_coarsening_and_refinement();
    }
}


void
distribute_dofs(DoFHandler<2> &dof_handler)
{
  static const FE_Q<2> finite_element(
    1); // degree of poly 1: bilinear element, 2:bi-quadratic element, ...
  // std::cout << "Value at (1,1) "
  //           << finite_element.shape_value(2, Point<2>(1, 1)) << "\n";


  dof_handler.distribute_dofs(
    finite_element); // we have associated a degree of freedom with a global
                     // number to each vertex

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());



  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);
  std::cout << "Number of entries per row in the sparsity pattern"
            << "\n";
  // SparseMatrix<double> system_matrix;
  // system_matrix.reinit(sparsity_pattern);

  for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
      std::cout << "Length of " << i << "-th row"
                << "\t" << sparsity_pattern.row_length(i) << "\n";
    }


  show_statistics(sparsity_pattern);
  std::ofstream out("sparsity_pattern1.svg");
  sparsity_pattern.print_svg(out);


  std::cout
    << "Printing all the non-zero entries for row 42 before renumbering: "
    << "\n";

  auto start_row_42 = sparsity_pattern.begin(42);
  auto end_row_42   = sparsity_pattern.end(42);
  for (auto it = start_row_42; it != end_row_42; ++it)
    {
      std::cout << (*it).column() << "\n";
    }
}



void
renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);


  std::cout << "Row length for Line 42: " << sparsity_pattern.row_length(42)
            << "\n";

  show_statistics(sparsity_pattern); // Bonus
  std::ofstream out("sparsity_pattern2.svg");
  sparsity_pattern.print_svg(out);

  std::cout
    << "Printing all the non-zero entries for row 42 after renumbering: "
    << "\n";
  for (auto it = sparsity_pattern.begin(42); it != sparsity_pattern.end(42);
       ++it)
    {
      std::cout << (*it).column() << "\n";
    }
}


void
make_grid_square(Triangulation<2> &square)
{
  GridGenerator::hyper_cube(square, -1, 1);
  square.refine_global(3);
}



int
main()
{
  Triangulation<2> triangulation;
  make_grid(triangulation);

  DoFHandler<2> dof_handler(triangulation);

  distribute_dofs(dof_handler);
  renumber_dofs(dof_handler);


  // //Test with a square
  // Triangulation<2> square;
  // make_grid_square(square);

  // DoFHandler<2> dof_handler_square(square);

  // distribute_dofs(dof_handler_square);
  // renumber_dofs(dof_handler_square);
}
