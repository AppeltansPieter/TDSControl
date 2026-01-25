#include "../src/utils/grid_iterator.hpp"

#include "../src/utils/ipow.hpp"

#include <gtest/gtest.h>

using namespace tds::utils;

TEST(grid_iterator, omega_iterator_one_dimensional)
{
    const ThetaRange range(1, 10);
    auto iter = range.begin();
    EXPECT_EQ(*iter, std::vector{0UL});
    ++iter;
    EXPECT_EQ(*iter, std::vector{1UL});
    ++iter;
    ++iter;
    ++iter;
    ++iter;
    EXPECT_EQ(*iter, std::vector{5UL});
    EXPECT_NE(iter, range.end());
    ++iter;
    EXPECT_EQ(iter, range.end());
}

TEST(grid_iterator, omega_iterator_two_dimensional)
{
    std::size_t count = 0;
    const ThetaRange range(2, 20);
    auto iter = range.begin();
    EXPECT_EQ(*iter, (std::vector{0UL, 0UL}));
    ++iter;
    EXPECT_EQ(*iter, (std::vector{1UL, 0UL}));
    count++;
    for (; count < 10; ++count) {
        ++iter;
    }
    EXPECT_EQ(*iter, (std::vector{10UL, 0UL}));
    ++iter;
    EXPECT_EQ(*iter, (std::vector{0UL, 1UL}));
    ++count;
    for (; iter != ThetaRange::end(); ++iter) {
        ++count;
    }
    EXPECT_EQ(iter, range.end());
    EXPECT_EQ(count, 20 * 11);
}

TEST(grid_iterator, omega_iterator_five_dimensional)
{
    std::size_t count = 0;
    const ThetaRange range(5, 6);
    auto iter = range.begin();
    EXPECT_EQ(*iter, (std::vector{0UL, 0UL, 0UL, 0UL, 0UL}));
    ++iter;
    EXPECT_EQ(*iter, (std::vector{1UL, 0UL, 0UL, 0UL, 0UL}));
    count++;
    for (; count < 3; ++count) {
        ++iter;
    }
    EXPECT_EQ(*iter, (std::vector{3UL, 0UL, 0UL, 0UL, 0UL}));
    ++iter;
    EXPECT_EQ(*iter, (std::vector{0UL, 1UL, 0UL, 0UL, 0UL}));
    ++count;
    for (; iter != ThetaRange::end(); ++iter) {
        ++count;
    }
    EXPECT_EQ(iter, range.end());
    EXPECT_EQ(count, ipow(6, 4) * 4);
}
