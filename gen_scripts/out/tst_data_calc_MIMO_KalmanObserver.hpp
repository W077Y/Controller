#pragma once
#ifndef TST_DATA_CALC_MIMO_KALMANOBSERVER_HPP_INCLUDED
#define TST_DATA_CALC_MIMO_KALMANOBSERVER_HPP_INCLUDED

#include <exmath.hpp>

#include <controller.hpp>

namespace test_data::parameter
{
  using parameter_t = controller::calculator::MIMO_KalmanObserver<float, 5, 2, 2>::parameter_t;
  extern parameter_t const param;
  extern exmath::matrix_t<float, 2, 2000> const y_meas;
  extern exmath::matrix_t<float, 2, 2000> const y;
  extern exmath::matrix_t<float, 2, 2000> const u;
  extern exmath::matrix_t<float, 5, 2000> const x;
} // namespace test_data::parameter

#endif

