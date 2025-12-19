#include <tdscontrol/tds.hpp>
#include "assertions.hpp"

namespace tds
{
    tds::tds(std::vector<double> A, std::vector<double> hA) 
    : m_A(std::move(A)), m_hA(std::move(hA))
    {
        TDS_CONTROL_PRECONDITION(m_A.size() == m_hA.size(), "Number of elements in A and hA do not match!");
    }
    #ifdef JULIA_BINDING
    tds::tds(jlcxx::ArrayRef<double> A, jlcxx::ArrayRef<double> hA)
    :m_A(A.begin(), A.end()), m_hA(hA.begin(), hA.end()) 
    { 
        TDS_CONTROL_PRECONDITION(m_A.size() == m_hA.size(), "Number of elements in A and hA do not match!");
    }
    #endif

    const std::vector<double> &tds::A() const {
        return m_A;
    }
    
    const std::vector<double> &tds::hA() const {
        return m_hA;
    }
    
} // namespace tds
