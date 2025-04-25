#pragma once
#ifndef CONTROLLER_HPP_INCLUDED
#define CONTROLLER_HPP_INCLUDED

#include <controller_Interface.hpp>
#include <exmath.hpp>

namespace controller::calculator
{
  template <typename T, std::size_t N> class SISO_KalmanIntController
  {
    using value_type                = std::remove_cvref_t<T>;
    static constexpr std::size_t nn = N;
    static constexpr std::size_t no = 1;
    static constexpr std::size_t ni = 1;
    static constexpr std::size_t nr = 1;

  public:
    struct parameter_t
    {
      static constexpr exmath::matrix_t<value_type, nn, nn> eye = exmath::matrix_t<value_type, nn, nn>::unit();

      exmath::matrix_t<value_type, nn, nn> const A_kal;
      exmath::matrix_t<value_type, nn, ni> const B_kal;
      exmath::matrix_t<value_type, no, nn> const C_kal;
      exmath::matrix_t<value_type, nn, nn> const Q_kal;
      exmath::matrix_t<value_type, no, no> const R_kal;
      exmath::matrix_t<value_type, ni, nn> const hT_kal;

      exmath::matrix_t<value_type, nr, 2 * ni + nr> const B_int;
      exmath::matrix_t<value_type, nr, nn> const          C_int;
      exmath::matrix_t<value_type, ni, nr> const          hT_int;
    };

    struct actuator_values_t
    {
      value_type u_star;
      value_type u_dash;
      value_type u;
    };

    SISO_KalmanIntController(parameter_t const& param)
        : m_param{ param }
    {
    }

    exmath::matrix_t<value_type, no, 1> calculate_y_kal() const { return this->m_param.C_kal * this->m_x_kal; }
    exmath::matrix_t<value_type, nr, 1> calculate_y_int() const { return this->m_param.C_int * this->m_x_kal; }
    exmath::matrix_t<value_type, ni, 1> calculate_u_star() const { return 0.0f - this->m_param.hT_kal * this->m_x_kal - this->m_param.hT_int * this->m_x_int; }
    exmath::matrix_t<value_type, nn, 1> const& get_x() const { return this->m_x_kal; }

    void correction_step(exmath::matrix_t<value_type, no, 1> const& y_meas)
    {
      auto const K = (this->m_pp * exmath::transpose(this->m_param.C_kal)) /
                     (this->m_param.C_kal * this->m_pp * exmath::transpose(this->m_param.C_kal) + this->m_param.R_kal);
      auto const p_ = (this->m_param.eye - K * this->m_param.C_kal) * this->m_pp;
      this->m_pp    = this->m_param.A_kal * p_ * exmath::transpose(this->m_param.A_kal) + this->m_param.Q_kal;

      this->m_x_kal += K * (y_meas - this->calculate_y_kal());
    }

    void update_step(value_type const& e, actuator_values_t const& input)
    {
      this->m_x_kal = this->m_param.A_kal * this->m_x_kal + this->m_param.B_kal * input.u;
      this->m_x_int += this->m_param.B_int * exmath::matrix_t<value_type, 2 * ni + nr, 1>{ input.u_dash, input.u_star, e };
    }

    void reset()
    {
      this->m_pp    = this->m_param.eye;
      this->m_x_kal = {};
      this->m_x_int = {};
    }

  private:
    parameter_t const& m_param;

    exmath::matrix_t<value_type, nn, nn> m_pp    = this->m_param.eye;
    exmath::matrix_t<value_type, nn, 1>  m_x_kal = {};
    exmath::matrix_t<value_type, nr, 1>  m_x_int = {};
  };
}    // namespace controller::calculator

#endif
