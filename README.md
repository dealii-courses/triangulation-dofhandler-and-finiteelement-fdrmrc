#  Lab 02 - Triangulation, DoFHandler, FiniteElement
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

1.  Setup

    - You will need the following programs:
      - Visual Studio Code with the `Remote Development` plugin installed
      - Paraview (ver > 5.4 would be better)
    -   Check that the `.devcontainer` directory contains a configuration that works for you. When you open the root directory of the repository, it should prompt you to reopen in container. If not, try explictily with the command `F1 -> Reopen folder in container`
    - If opening in container worked, you should be able to run a terminal (`CTRL + ~`) and see something like
    ```
    To run a command as administrator (user "root"), use "sudo <command>".
    See "man sudo_root" for details.

    dealii@3f25621373bd:/workspaces/02_tria_dh_fem$ 
    ```
    - The folder contains a `CMakeLists.txt` file. VSCode should detect this, and provide commands to configure and build your application, at the bottom of the terminal window

2.  Tasks for lab-02:

    1.  Read the documentation of step-1 deal.II at\
        -   <https://www.dealii.org/current/doxygen/deal.II/step_1.html>
        -   <https://www.dealii.org/current/doxygen/deal.II/step_2.html>

    2.  Compile and run `step-1` inside the VSCode container and look at the output.

    3.  Comment out the ``.set_manifold(0, ...)`` line in
        ``second_grid()``. What happens now?

    4.  Create an image of an L-shape domain (add a function third_grid() to
        step-1) with one global refinement.

    5.  Now change the output format of the previous example to vtk and open
    	the new file in paraview.

    6.  Refine the L-shaped mesh adaptively around the re-entrant corner
        several times but with a twist: refine all cells with the distance
        between the center of the cell and re-entrant corner is smaller than
        1/3.

    7.  Output mesh two as an svg file instead of eps. Open it in a
        browser to display it (firefox for example).

    8.  Create a helper function that takes a reference to a
        Triangulation and returns a tuple with  the following information: number of levels, number of cells, number of active cells. Test this with all of your meshes.

    9.  Generate a circle using ``GridGenerator::hyper_ball()`` in 2d: use a
        SphericalManifold everywhere, only on the boundary, or on all cells
        except the center cell and refine the mesh globally twice.

    10. Go into ``second_grid()`` and add a `triangulation.reset_manifold(0)` right before the refinment loop. Explain in a comment of the code what is happening.

    11. Bonus: Create a mesh that represents the surface of a torus and refine
        it 2 times globally. Output to vtk format and check the output. Note
        that your Triangulation needs to be of type ``Triangulation<2,3>``,
        which we will discuss later this week.

    12. Bonus: Take a look at step-49 and read the included .msh file in your modified step-1 program.

    13.  Read the documentation of step-2 deal.II at\
    -   <https://www.dealii.org/current/doxygen/deal.II/step_2.html>

    14. Run step-2. Look at the sparsity patterns in firefox.

    15.  How does the pattern change if you increase the polynomial degree
         from 1 to 2 or to 3?

    16.  How does the pattern change if you use a globally refined (say 3
         times) unit square?

    16.  Are these patterns symmetric? Why/why not?

    17.  How many entries per row in the sparsity pattern do you expect for a
         Q1 element (assuming four cells are around each vertex)? Check that
         this is true for the mesh in b) (look for `row_length(i)` and output
         them for each row). Can you construct a 2d mesh (without hanging
         nodes) that has a row with more entries?

    18.  How many entries per row in the sparsity pattern are there for Q2 and
        Q3 elements, again assuming four cells around each vertex?

    19.  Print all entries for row 42 for the original renumbered sparsity
         pattern.

    20.  Bonus: Compute and output statistics like the number of
         unknowns, bandwidth of the sparsity pattern, average number of
         entries per row, and fill ratio.





