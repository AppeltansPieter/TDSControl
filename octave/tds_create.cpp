#include <octave/oct.h>
#include <octave/dMatrix.h>
#include <octave/dRowVector.h>
#include <octave/Cell.h>

#include "octave_TDS.hpp"


DEFUN_DLD (tds_create, args, nargout,
           "tds_create")
{
    octave_value_list retval;

    Cell A_arg = args(0).cell_value();
    RowVector hA_arg = args(1).row_vector_value();
    // TODO: check A_arg.numel() == hA_arg.numel();
    std::vector<double> A;
    A.reserve(A_arg.numel());
    std::vector<double> hA;
    hA.reserve(hA_arg.numel()); 
    for (std::size_t i = 0; i < hA_arg.numel(); i++) {
        A.push_back(A_arg(i).double_value());
        hA.push_back(hA_arg(i));
    }
    auto *res = new octave_TDS(A, hA);
    retval.append(octave_value(res));

    return retval;
}