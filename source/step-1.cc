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
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>

using namespace dealii;


std::tuple<unsigned int, unsigned int, unsigned int>
get_info_triangulation(const Triangulation<2> &tria)
{
  return std::make_tuple<unsigned int, unsigned int>(tria.n_active_cells(),
                                                     tria.n_cells(),
                                                     tria.n_levels());
}



void
first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  std::ofstream out("grid-1.vtk");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-1.vtk" << std::endl;

  auto info_tria = get_info_triangulation(triangulation);
  std::cout << std::get<0>(info_tria) << "\t" << std::get<1>(info_tria) << "\t"
            << std::get<2>(info_tria) << "\n";
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
      triangulation.reset_manifold(
        0); // Reset all parts of the triangulation, regardless of their
            // manifold_id, to use a FlatManifold object.
      triangulation.execute_coarsening_and_refinement();
    }


  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  auto info_tria = get_info_triangulation(triangulation);

  std::cout << "Grid written to grid-2.svg" << std::endl;
  std::cout << std::get<0>(info_tria) << "\t" << std::get<1>(info_tria) << "\t"
            << std::get<2>(info_tria) << "\n";
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
  auto info_tria = get_info_triangulation(triangulation);
  std::cout << std::get<0>(info_tria) << "\t" << std::get<1>(info_tria) << "\t"
            << std::get<2>(info_tria) << "\n";

  std::cout << "Now refine the grid globally" << std::endl;
  const unsigned int initial_global_refinement = 2;

  triangulation.refine_global(initial_global_refinement);
  std::ofstream out_refined("grid-4.svg");

  grid_out.write_svg(triangulation, out_refined);
  std::cout << "Grid written to grid-4.svg" << std::endl;


  std::cout << "Refine adaptively around the re-entrant corner"
            << "\n";
  const Point<2>     corner(0.0, 0.0);
  const unsigned int no_refinements{4};
  for (unsigned int step = 0; step < no_refinements; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_corner =
                corner.distance(cell->center(v));
              // std::cout << cell->center(v) << cell->center() << "\n"; //check
              // about syntax
              if (std::fabs(distance_from_corner) <= (1.0) / (3.0))
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }
      triangulation.execute_coarsening_and_refinement();
    }


  std::ofstream out_refined_locally("grid-5.svg");

  grid_out.write_svg(triangulation, out_refined_locally);
  std::cout << "Grid written to grid-5.svg" << std::endl;


  auto info_tria_L = get_info_triangulation(triangulation);
  std::cout << std::get<0>(info_tria_L) << "\t" << std::get<1>(info_tria_L)
            << "\t" << std::get<2>(info_tria_L) << "\n";
}



void
fourth_grid()
{
  Triangulation<2> triangulation;
  const double     radius = 1.0;
  const Point<2>   center(0.0, 0.0); // center of the ball


  const SphericalManifold<2> manifold(center);
  GridGenerator::hyper_ball(
    triangulation,
    center,
    radius,
    true); // spherical manifold set to true on the bdary
  std::ofstream out("grid-Spherical_Manifold_bdary.svg");
  GridOut       grid_out;

  std::cout << "Using Spherical Manifold on the boundary "
            << "\n";
  triangulation.reset_all_manifolds();
  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold(0, manifold);
  triangulation.refine_global(2);

  grid_out.write_svg(triangulation, out);


  //////////////EVERYWHERE
  Triangulation<2> triangulation_e;
  GridGenerator::hyper_ball(triangulation_e,
                            center,
                            radius,
                            true); // spherical manifold set to true

  std::ofstream out_e("grid-Spherical_Manifold_everywhere.svg");
  GridOut       grid_out_e;

  std::cout << "Using Spherical Manifold everywhere "
            << "\n";
  triangulation_e.reset_all_manifolds();
  // reenable the manifold:
  triangulation_e.set_all_manifold_ids(0);
  triangulation_e.set_manifold(0, manifold);
  triangulation_e.refine_global(2);

  grid_out_e.write_svg(triangulation_e, out_e);



  //////////EXCEPT THE CENTER
  Triangulation<2> triangulation_c;
  GridGenerator::hyper_ball(triangulation_c, center, radius, true);

  std::ofstream  out_c("grid-Spherical_Manifold_everywhere_but_center.svg");
  GridOut        grid_out_c;
  const Point<2> mesh_center;
  for (const auto &cell : triangulation_c.active_cell_iterators())
    {
      if (mesh_center.distance(cell->center()) > cell->diameter() / 10)
        {
          cell->set_all_manifold_ids(0);
        }
    }
  std::cout << "Using Spherical Manifold everywhere but the center"
            << "\n";
  triangulation_c.refine_global(4);

  grid_out_c.write_svg(triangulation_c, out_c);
}



int
main()
{
  first_grid();
  second_grid();
  third_grid();
  fourth_grid();
}
