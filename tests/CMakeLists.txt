  message(STATUS "Building tests")
  include(CTest)
  include(Sanitizers)
  include(Warnings)
  include(AddTest)

  enable_testing()
  find_package(GTest REQUIRED)


  enable_sanitizers(tdscontrol)
  enable_warnings(tdscontrol)
  target_compile_definitions(tdscontrol PUBLIC TDS_TESTING)

  add_tdscontrol_test(test_roots "test_roots.cpp")
  add_tdscontrol_test(test_newton "newtons_method.cpp")
  add_tdscontrol_test(test_sparse_spectral_discretization
                      "sparse_spectral_discretization.cpp")
  add_tdscontrol_test(test_grid_iterator "utils/grid_iterator.cpp")
  add_tdscontrol_test(test_splines "utils/cubic_spline.cpp")
  add_tdscontrol_test(test_compute_N_rhp "compute_N_rhp.cpp")
