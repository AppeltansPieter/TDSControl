#include "assertions.hpp"

#include <cmath>
#include <tdscontrol/tds.hpp>
#include <utility>
#include <vector>

namespace tds {
tds::tds(std::vector<Eigen::MatrixXd> A, std::vector<double> hA) : m_A(std::move(A)), m_hA(std::move(hA))
{
    TDS_CONTROL_PRECONDITION(m_A.size() == m_hA.size(), "Number of elements in A and hA do not match!");
    TDS_CONTROL_PRECONDITION(std::ranges::all_of(m_A, [](const Eigen::MatrixXd &m1) { return m1.rows() == m1.cols(); }),
                             "Non-square matrix in A.");
    TDS_CONTROL_PRECONDITION(
        std::ranges::adjacent_find(m_A, [](const Eigen::MatrixXd &m1,
                                           const Eigen::MatrixXd &m2) { return m1.rows() != m2.rows(); }) == m_A.end(),
        "Dimensions of the matrices in A do not match!");
    // Test
}

} // namespace tds
