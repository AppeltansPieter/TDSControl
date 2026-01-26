#pragma once
#include <Eigen/Core>
#include <vector>

namespace tds {

// TODO: deal with commensurate delays
class tds {
  public:
    tds(std::vector<Eigen::MatrixXd> A, std::vector<double> hA);

    const std::vector<Eigen::MatrixXd> &A() const
    {
        return m_A;
    }
    const Eigen::MatrixXd &A(const std::size_t k) const
    {
        return m_A.at(k);
    }
    Eigen::MatrixXd &A(const std::size_t k)
    {
        return m_A.at(k);
    }

    const std::vector<double> &hA() const
    {
        return m_hA;
    }

    const double &hA(const std::size_t k) const
    {
        return m_hA.at(k);
    }

    double &hA(const std::size_t k)
    {
        return m_hA.at(k);
    }

    std::size_t mA() const
    {
        return m_A.size();
    }

    Eigen::Index n() const
    {
        return m_A.empty() ? 0 : m_A.at(0).rows();
    }

  private:
    std::vector<Eigen::MatrixXd> m_A;
    std::vector<double> m_hA;
};

} // namespace tds
