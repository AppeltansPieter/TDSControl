#pragma once

#include <cstdlib>
#include <vector>

namespace tds::utils {

/**
 * @brief Range object to iterate over the grid [0, -(p -1) / p *pi, ..., -(p-1) / p * p] -> [pi, pi, ..., pi]
 */
class ThetaRange : public std::ranges::view_interface<ThetaRange> {
  public:
    class iterator {
      public:
        using value_type = std::vector<std::size_t>;
        using reference = const value_type &;
        using pointer = const value_type *;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag;

        iterator(const std::size_t m, const std::size_t p, const bool end = false) : omega_(m), p_(p), end_(end) {}
        explicit iterator(const bool end) : omega_(0), p_(0), end_(end) {}

        reference operator*() const
        {
            return omega_;
        }
        pointer operator->() const
        {
            return &omega_;
        }

        iterator &operator++();
        bool operator==(const iterator &other) const
        {
            return end_ == other.end_ && (end_ || omega_ == other.omega_);
        }
        bool operator!=(const iterator &other) const
        {
            return !(*this == other);
        }

      private:
        // TODO: use union to save space?
        std::vector<std::size_t> omega_;
        std::size_t p_;
        bool end_; // TODO: Remove?
    };
    ThetaRange(const std::size_t m, const std::size_t p) : m_(m), p_(p) {}
    iterator begin() const
    {
        return {m_, p_};
    }
    static iterator end()
    {
        return iterator{true};
    }

  private:
    std::size_t m_;
    std::size_t p_;
};

} // namespace tds::utils