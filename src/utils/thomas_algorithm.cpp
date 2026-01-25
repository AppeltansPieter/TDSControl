#include "thomas_algorithm.hpp"

#include <cassert>

void tds::utils::thomas_algorithm(const std::vector<double> &sub, std::vector<double> &diag,
                                  const std::vector<double> &super, std::vector<double> &rhs)
{
    assert(sub.size() == diag.size() - 1);
    assert(super.size() == diag.size() - 1);
    assert(rhs.size() == diag.size());

    for (size_t i = 1; i < diag.size(); ++i) {
        const double w = sub[i - 1] / diag[i - 1];
        diag[i] -= w * super[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }
    rhs.back() = rhs.back() / diag.back();

    for (size_t i = 1; i < rhs.size(); ++i) {
        const size_t ind = rhs.size() - 1 - i;
        rhs[ind] = (rhs[ind] - super[ind] * rhs[ind + 1]) / diag[ind];
    }
}