#include <iomanip>
#include <iostream>
#include <tdscontrol/roots.hpp>
#include <tdscontrol/tds.hpp>

int main(int argc, char const *argv[])
{
    tds::tds sys({1, -1}, {0, 1});
    const auto roots = tds::roots(sys, 15);
    std::cout << "Eigenvalues: ";
    std::cout << std::fixed << std::setprecision(16);
    for (const auto &root : roots) {
        std::cout << root << ", ";
    }
    std::cout << "\n";

    return 0;
}
