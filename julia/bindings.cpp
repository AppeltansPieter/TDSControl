#include "jlcxx/jlcxx.hpp"
#include <jlcxx/stl.hpp>     // Required for std::vector
#include <tdscontrol/tds.hpp>
#include <tdscontrol/roots.hpp>

using namespace std::literals::complex_literals;


JLCXX_MODULE define_julia_module(jlcxx::Module &mod){
    mod.add_type<tds::tds>("TDS")
        .constructor<jlcxx::ArrayRef<double>, jlcxx::ArrayRef<double>>()
        .method("get_mA", &tds::tds::mA)
        .method("get_n", &tds::tds::n);

    mod.method("roots",[](const tds::tds &sys, std::size_t N) { 
        auto tmp = tds::roots(sys, N);
        jlcxx::Array<std::complex<double>> res{};
        for (auto &elem: tmp) {
            res.push_back(elem);
        }
        return res;
    }
);
}