#pragma once
#include <vector>

namespace tds::utils {

/**
 * @brief Use the Thomas algorithm to solve a tridiagonal system.
 *
 * @param sub [in] Vector containing the subdiagonal elements
 * @param diag [in] Vector containing the diagonal elements, overwritten in the computation process.
 * @param super [in] Vector containing the superdiagonal elements
 * @param rhs [in/out] Vector containing the right hand-side elements, overwritten with the solution.
 */
void thomas_algorithm(const std::vector<double> &sub, std::vector<double> &diag, const std::vector<double> &super,
                      std::vector<double> &rhs);
} // namespace tds::utils