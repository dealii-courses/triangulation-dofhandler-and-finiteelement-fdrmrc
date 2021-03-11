#include <deal.II/base/point.h>

#include <gtest/gtest.h>

using namespace dealii;

// Declare functions in step-1 and step-2
void
first_grid();
void
second_grid();


TEST(Step1, first_grid)
{
  first_grid();
}

TEST(Step1, second_grid)
{
  second_grid();
}


int
main(int argc, char *argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
