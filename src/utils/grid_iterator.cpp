#include "grid_iterator.hpp"

#include <filesystem>

namespace tds::utils {

ThetaRange::iterator &ThetaRange::iterator::operator++()
{
    if (end_) {
        return *this;
    }
    omega_[0]++;
    std::size_t l = 0;
    while (l < omega_.size() - 1 && omega_[l] > (l > 0 ? p_ - 1 : p_ / 2)) {
        omega_[l] = 0;
        l++;
        omega_[l] += 1;
    }
    if (omega_.back() > ((omega_.size() == 1) ? p_ / 2 : p_ - 1)) {
        end_ = true;
    }
    return *this;
}
} // namespace tds::utils
