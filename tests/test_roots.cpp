#include <gtest/gtest.h>
#include <tdscontrol/roots.hpp>
#include <tdscontrol/tds.hpp>

GTEST_TEST(roots, scalar)
{
    tds::tds sys({1., -1.}, {0., 1.});
    auto roots = tds::roots(sys, 10);
    EXPECT_EQ(roots.size(), 11);
}
