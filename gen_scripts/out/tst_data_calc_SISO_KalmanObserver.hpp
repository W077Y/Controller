#pragma once
#ifndef TST_DATA_CALC_SISO_KALMANOBSERVER_HPP_INCLUDED
#define TST_DATA_CALC_SISO_KALMANOBSERVER_HPP_INCLUDED

#include <exmath.hpp>

#include <controller.hpp>

namespace test_data::parameter
{
  using parameter_t = controller::calculator::SISO_KalmanObserver<float, 3>::parameter_t;
  extern parameter_t const param;
  extern exmath::matrix_t<float, 1, 5000> const y_meas;
  extern exmath::matrix_t<float, 1, 5000> const y;
  extern exmath::matrix_t<float, 1, 5000> const u;
  extern exmath::matrix_t<float, 3, 5000> const x;
} // namespace test_data::parameter

#endif

