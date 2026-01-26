#include "../../src/utils/cubic_spline.hpp"

#include <gtest/gtest.h>

using namespace tds::utils;

TEST(Splines, cubic_spline)
{
    const CubicSpline s{0., 1., std::array{2., 3., 5., 6., 8., 4., 1., -5.}};
    EXPECT_FLOAT_EQ(s(0.), 2.);
    EXPECT_FLOAT_EQ(s(3.), 6.);
    EXPECT_FLOAT_EQ(s(7.), -5.);
    EXPECT_FLOAT_EQ(s(0.5), 2.352885606320852);
    EXPECT_FLOAT_EQ(s(3.5), 7.406690140845071);
    EXPECT_FLOAT_EQ(s(6.5), -1.625772930264514);
    EXPECT_FLOAT_EQ(s(4.3), 7.232213328753007);
    EXPECT_FLOAT_EQ(s(2.7), 5.539603572655445);
}