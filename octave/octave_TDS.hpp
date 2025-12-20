#pragma once

#include <octave/oct.h>
#include <octave/dRowVector.h>
#include <octave/Cell.h>

#include <memory>

#include <tdscontrol/tds.hpp>


class octave_TDS: public octave_base_value {
    public:
        octave_TDS() = default;

        octave_TDS(const std::vector<double> &A, const std::vector<double> &hA) : m_tds(std::make_shared<tds::tds>(A, hA)) {
        }

        octave_TDS(const octave_TDS &other)
        : m_tds(std::make_shared<tds::tds>(*other))
        {
        }

        const octave_TDS &operator=(const octave_TDS &other) {
            m_tds = std::make_shared<tds::tds>(*other);
            return *this;
        }

        bool is_defined() const override { return m_tds.operator bool(); }

        octave_base_value *clone() const override { 
            return new octave_TDS(*this); 
        }

        octave_base_value *empty_clone() const override { 
            return new octave_TDS(); }

        const tds::tds &operator*() const {
            return *m_tds;
        }

        tds::tds &operator*() {
            return *m_tds;
        }

        void print(std::ostream& os, bool) override {
            os << "TDS(" << m_tds.get() << ")";
        }

        
    private:
        std::shared_ptr<tds::tds> m_tds;
        DECLARE_OV_BASE_TYPEID_FUNCTIONS_AND_DATA 
};

DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (octave_TDS,
                                     "tds",
                                     "octave_TDS");