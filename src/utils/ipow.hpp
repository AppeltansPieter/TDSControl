#pragma once
#include <cassert>
#include <cstdlib>
#include <limits>
#include <type_traits>

/**
 * @brief Compute the absolute value of an integer value.
 *
 * @tparam T  An integer type
 * @param x The value to take the absolute value from
 * @return The absolute value
 */
template <typename T> constexpr T iabs(T x)
{
    if constexpr (std::is_signed_v<T>) {
        return std::abs(x);
    }
    else {
        return x;
    }
}

/**
 * @brief Raise an integer to an integer power.
 *
 * @note The function does not account for integer overflow.
 *
 * @tparam T An integer type
 * @param base The number to raise to the power
 * @param exp
 * @return
 */
template <typename T> constexpr T ipow(T base, unsigned exp)
{
    T result = 1;
    while (exp > 0) {
        if (exp & 1) {
            assert(iabs(result) < std::numeric_limits<T>::max() / iabs(base));
            result *= base;
        }
        assert(iabs(base) < std::numeric_limits<T>::max() / iabs(base));
        base *= base;
        exp >>= 1;
    }
    return result;
}
