#pragma once
#include <Eigen/Dense>
#include <functional>
#include <iostream>

namespace tds::internal {

struct Tolerance {
    double value;
};
struct MaxIter {
    std::size_t value;
};

class NewtonsMethod {
  public:
    explicit NewtonsMethod(const Eigen::Index n) : m_fun(n), m_step(n), m_jacobian(n, n), m_qr(m_jacobian) {}
    void operator()(Eigen::VectorXcd &x,
                    const std::function<void(const Eigen::VectorXcd &, Eigen::VectorXcd &, Eigen::MatrixXcd &)> &f,
                    const Tolerance tol = Tolerance(1e-10), const MaxIter max_iter = MaxIter(10))
    {
#ifdef TDS_TESTING
        m_residuals.clear();
#endif
        std::size_t iter = 0;
        // TODO: Compute Jacobian in separate function?
        f(x, m_fun, m_jacobian);
        double residual = m_fun.norm();
        while (residual > tol.value) {
#ifdef TDS_TESTING
            m_residuals.push_back(residual);
#endif
            m_qr.compute(m_jacobian);
            m_step = m_qr.solve(m_fun);
            x.noalias() -= m_step;
            double relative_step_size = m_step.norm() / x.norm();
            if (relative_step_size < tol.value) {
                break;
            }
            if (++iter >= max_iter.value) {
                break;
            }
            f(x, m_fun, m_jacobian);
            residual = m_fun.norm();
        }
#ifdef TDS_TESTING
        m_residuals.push_back(residual);
#endif
    }
#ifdef TDS_TESTING
    const std::vector<double> &get_residuals() const
    {
        return m_residuals;
    }
#endif
  private:
    Eigen::VectorXcd m_fun;
    Eigen::VectorXcd m_step;
    Eigen::MatrixXcd m_jacobian;
    Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::MatrixXcd>> m_qr;
#ifdef TDS_TESTING
    std::vector<double> m_residuals;
#endif
};
} // namespace tds::internal
