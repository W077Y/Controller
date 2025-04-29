#pragma once
#ifndef CONTROLLER_CALCULATORS_HPP_INCLUDED
#define CONTROLLER_CALCULATORS_HPP_INCLUDED

#include <exmath.hpp>

namespace controller::calculator
{
  template <typename T, std::size_t N> class SISO_KalmanObserver
  {
    using value_type                = std::remove_cvref_t<T>;
    static constexpr std::size_t nn = N;
    static constexpr std::size_t no = 1;
    static constexpr std::size_t ni = 1;

  public:
    struct parameter_t
    {
      static constexpr exmath::matrix_t<value_type, nn, nn> eye = exmath::matrix_t<value_type, nn, nn>::unit();

      exmath::matrix_t<value_type, nn, nn> const A;
      exmath::matrix_t<value_type, nn, ni> const B;
      exmath::matrix_t<value_type, no, nn> const C;
      exmath::matrix_t<value_type, nn, nn> const Q;
      exmath::matrix_t<value_type, no, no> const R;
    };

    SISO_KalmanObserver(parameter_t const& param)
        : m_param{ param }
    {
    }

    exmath::matrix_t<value_type, no, 1>        get_y() const { return this->m_param.C * this->m_x; }
    exmath::matrix_t<value_type, nn, 1> const& get_x() const { return this->m_x; }

    void correction_step(exmath::matrix_t<value_type, no, 1> const& y_meas)
    {
      auto const K  = (this->m_pp * exmath::transpose(this->m_param.C)) / (this->m_param.C * this->m_pp * exmath::transpose(this->m_param.C) + this->m_param.R);
      auto const p_ = (this->m_param.eye - K * this->m_param.C) * this->m_pp;
      this->m_pp    = this->m_param.A * p_ * exmath::transpose(this->m_param.A) + this->m_param.Q;

      this->m_x += K * (y_meas - this->get_y());
    }

    void update_step(exmath::matrix_t<value_type, ni, 1> const& input) { this->m_x = this->m_param.A * this->m_x + this->m_param.B * input; }

    void reset()
    {
      this->m_pp = this->m_param.eye;
      this->m_x  = {};
    }

  private:
    parameter_t const& m_param;

    exmath::matrix_t<value_type, nn, nn> m_pp = this->m_param.eye;
    exmath::matrix_t<value_type, nn, 1>  m_x  = {};
  };

  template <typename T, std::size_t N, std::size_t I, std::size_t O> class MIMO_KalmanObserver
  {
    using value_type                = std::remove_cvref_t<T>;
    static constexpr std::size_t nn = N;
    static constexpr std::size_t ni = I;
    static constexpr std::size_t no = O;

  public:
    struct parameter_t
    {
      static constexpr exmath::matrix_t<value_type, nn, nn> eye = exmath::matrix_t<value_type, nn, nn>::unit();

      exmath::matrix_t<value_type, nn, nn> const A;
      exmath::matrix_t<value_type, nn, ni> const B;
      exmath::matrix_t<value_type, no, nn> const C;
      exmath::matrix_t<value_type, nn, nn> const Q;
      exmath::matrix_t<value_type, no, no> const R;
    };

    MIMO_KalmanObserver(parameter_t const& param)
        : m_param{ param }
    {
    }

    exmath::matrix_t<value_type, no, 1>        get_y() const { return this->m_param.C * this->m_x; }
    exmath::matrix_t<value_type, nn, 1> const& get_x() const { return this->m_x; }

    void correction_step(exmath::matrix_t<value_type, no, 1> const& y_meas)
    {
      auto const K  = (this->m_pp * exmath::transpose(this->m_param.C)) * exmath::inv(this->m_param.C * this->m_pp * exmath::transpose(this->m_param.C) + this->m_param.R);
      auto const p_ = (this->m_param.eye - K * this->m_param.C) * this->m_pp;
      this->m_pp    = this->m_param.A * p_ * exmath::transpose(this->m_param.A) + this->m_param.Q;

      this->m_x += K * (y_meas - this->get_y());
    }

    void update_step(exmath::matrix_t<value_type, ni, 1> const& input) { this->m_x = this->m_param.A * this->m_x + this->m_param.B * input; }

    void reset()
    {
      this->m_pp = this->m_param.eye;
      this->m_x  = {};
    }

  private:
    parameter_t const& m_param;

    exmath::matrix_t<value_type, nn, nn> m_pp = this->m_param.eye;
    exmath::matrix_t<value_type, nn, 1>  m_x  = {};
  };

  template <typename T, std::size_t N> class SISO_KalmanIntController
  {
    using value_type                = std::remove_cvref_t<T>;
    static constexpr std::size_t nn = N;
    static constexpr std::size_t no = 1;
    static constexpr std::size_t ni = 1;
    static constexpr std::size_t nr = 1;

  public:
    using observer_t = SISO_KalmanObserver<value_type, N>;

    struct parameter_t
    {
      observer_t::parameter_t const              observer;
      exmath::matrix_t<value_type, ni, nn> const hT_kal;

      exmath::matrix_t<value_type, nr, 2 * ni + nr> const B_int;
      exmath::matrix_t<value_type, nr, nn> const          C_int;
      exmath::matrix_t<value_type, ni, nr> const          hT_int;
    };

    SISO_KalmanIntController(parameter_t const& param)
        : m_param{ param }
    {
    }

    exmath::matrix_t<value_type, no, 1> calculate_y_kal() const { return this->m_kal.get_y(); }
    exmath::matrix_t<value_type, nr, 1> calculate_y_int() const { return this->m_param.C_int * this->m_kal.get_x(); }
    exmath::matrix_t<value_type, ni, 1> calculate_u_star() const
    {
      return 0.0f - this->m_param.hT_kal * this->m_kal.get_x() - this->m_param.hT_int * this->m_x_int;
    }
    exmath::matrix_t<value_type, nn, 1> const& get_x() const { return this->m_kal.get_x(); }

    void correction_step(exmath::matrix_t<value_type, no, 1> const& y_meas) { return this->m_kal.correction_step(y_meas); }

    void update_observer(exmath::matrix_t<value_type, ni, 1> const& u) { return this->m_kal.update_step(u); }

    void update_integrator(value_type const& e, value_type const& u_star, value_type const& u_dash)
    {
      this->m_x_int += this->m_param.B_int * exmath::matrix_t<value_type, 2 * ni + nr, 1>{ u_dash, u_star, e };
    }

    void reset_observer() { this->m_kal.reset(); }

    void reset_integrator() { this->m_x_int = {}; }

  private:
    parameter_t const& m_param;

    observer_t                          m_kal{ this->m_param.observer };
    exmath::matrix_t<value_type, nr, 1> m_x_int = {};
  };

}    // namespace controller::calculator

#endif
