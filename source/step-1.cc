/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * based on deal.II step-1
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>

using namespace dealii;


void
first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  std::ofstream out("grid-1.vtu");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.vtu" << std::endl;
}



void
second_grid()
{
  Triangulation<2> triangulation;

  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
  for (unsigned int step = 0; step < 5; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }


  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.svg" << std::endl;
}



void
third_grid()
{
  Triangulation<2> triangulation;
  const double     left  = -1.0;
  const double     right = +1.0;
  GridGenerator::hyper_L(triangulation, left, right, true);
  std::ofstream out("grid-3.svg");
  GridOut       grid_out;

  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-3.svg" << std::endl;

  std::cout << "Now refine the grid globally" << std::endl;
  const unsigned int initial_global_refinement = 2;

  triangulation.refine_global(initial_global_refinement);
  std::ofstream out_refined("grid-4.svg");

  grid_out.write_svg(triangulation, out_refined);
  std::cout << "Grid written to grid-4.svg" << std::endl;


  std::cout << "Refine adaptively around the re-entrant corner"
            << "\n";
  const Point<2>     corner(0.0, 0.0);
  const unsigned int no_refinements{6};
  for (unsigned int step = 0; step < no_refinements; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_corner =
                corner.distance(cell->center(v));

              if (std::fabs(distance_from_corner) <= (1.0) / (3.0))
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }
    }
  triangulation.execute_coarsening_and_refinement();


  std::ofstream out_refined_locally("grid-5.svg");

  grid_out.write_svg(triangulation, out_refined_locally);
  std::cout << "Grid written to grid-5.svg" << std::endl;
}



int
main()
{
  first_grid();
  second_grid();
  third_grid();
}
