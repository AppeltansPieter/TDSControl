#include <octave/oct.h>
#include <octave/CMatrix.h>

#include "octave_TDS.hpp"

#include <tdscontrol/roots.hpp>

DEFUN_DLD (tds_roots, args, nargout,
           "Compute the characteristic roots of a time-delay system.")
{
    int nargin = args.length ();
    octave_value_list retval;
    if (nargin != 2) {
        print_usage();
    } else {
        octave_base_value *obj = args(0).internal_rep();
        octave_TDS *tds_wrapper = dynamic_cast<octave_TDS *>(obj);
        if (!tds_wrapper){
          error("Expected octave_TDS object.");
        }
        if (!(*tds_wrapper).is_defined())
          error("Expected octave_TDS object with non void pointer.");
        tds::tds &sys = *(*tds_wrapper);
        std::size_t N = args(1).int_value();
        std::vector<std::complex<double>> res = tds::roots(sys, N);
        ComplexMatrix r(res.size(), 1);
        std::copy(res.cbegin(), res.cend(), r.fortran_vec());
        return retval.append(r);
    }

  return retval; 
}