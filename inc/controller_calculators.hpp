#pragma once
#ifndef CONTROLLER_CALCULATORS_HPP_INCLUDED
#define CONTROLLER_CALCULATORS_HPP_INCLUDED

#include <exmath.hpp>

namespace controller::calculator
{
  template <typename T, std::size_t N, std::size_t I, std::size_t O> class MIMO_System
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
      exmath::matrix_t<value_type, nn, 1> const  X0;
    };

    MIMO_System(parameter_t const& param)
        : m_param{ param }
    {
    }

    exmath::matrix_t<value_type, no, 1>        get_y() const { return this->m_param.C * this->m_x; }
    exmath::matrix_t<value_type, nn, 1> const& get_x() const { return this->m_x; }
    exmath::matrix_t<value_type, nn, 1>&       get_x() { return this->m_x; }

    void update_step(exmath::matrix_t<value_type, ni, 1> const& input) { this->m_x = this->m_param.A * this->m_x + this->m_param.B * input; }

    void reset() { this->m_x = this->m_param.X0; }

  private:
    parameter_t const& m_param;

    exmath::matrix_t<value_type, nn, 1> m_x = {};
  };

  template <typename T, std::size_t N> using SISO_System = MIMO_System<T, N, 1, 1>;

  template <typename T, std::size_t N, std::size_t I, std::size_t O> class MIMO_KalmanObserver
  {
    using value_type                = std::remove_cvref_t<T>;
    static constexpr std::size_t nn = N;
    static constexpr std::size_t ni = I;
    static constexpr std::size_t no = O;

  public:
    using system_t = MIMO_System<value_type, N, I, O>;

    struct parameter_t
    {
      static constexpr exmath::matrix_t<value_type, nn, nn> eye = exmath::matrix_t<value_type, nn, nn>::unit();

      system_t::parameter_t const                sys;
      exmath::matrix_t<value_type, nn, nn> const Q;
      exmath::matrix_t<value_type, no, no> const R;
      exmath::matrix_t<value_type, nn, nn> const P0;
    };

    MIMO_KalmanObserver(parameter_t const& param)
        : m_param{ param }
    {
    }

    exmath::matrix_t<value_type, no, 1>        get_y() const { return this->m_sys.get_y(); }
    exmath::matrix_t<value_type, nn, 1> const& get_x() const { return this->m_sys.get_x(); }

    void correction_step(exmath::matrix_t<value_type, no, 1> const& y_meas)
    {
      auto const K  = this->calculate_k();
      auto const p_ = (this->m_param.eye - K * this->m_param.sys.C) * this->m_pp;
      this->m_pp    = this->m_param.sys.A * p_ * exmath::transpose(this->m_param.sys.A) + this->m_param.Q;

      this->m_sys.get_x() += K * (y_meas - this->get_y());
    }

    void update_step(exmath::matrix_t<value_type, ni, 1> const& input) { return this->m_sys.update_step(input); }

    void reset()
    {
      this->m_sys.reset();
      this->m_pp = this->m_param.P0;
    }

  private:
    exmath::matrix_t<value_type, nn, no> calculate_k() const
    {
      if constexpr (no == 1)
      {
        return (this->m_pp * exmath::transpose(this->m_param.sys.C)) /
               (this->m_param.sys.C * this->m_pp * exmath::transpose(this->m_param.sys.C) + this->m_param.R);
      }
      else
      {
        return (this->m_pp * exmath::transpose(this->m_param.sys.C)) *
               exmath::inv(this->m_param.sys.C * this->m_pp * exmath::transpose(this->m_param.sys.C) + this->m_param.R);
      }
    }

    parameter_t const&                   m_param;
    system_t                             m_sys = { this->m_param.sys };
    exmath::matrix_t<value_type, nn, nn> m_pp  = { this->m_param.P0 };
  };

  template <typename T, std::size_t N> using SISO_KalmanObserver = MIMO_KalmanObserver<T, N, 1, 1>;

  template <typename T, std::size_t N, std::size_t I, std::size_t O, std::size_t R> class MIMO_KalmanIntegralController
  {
    using value_type                = std::remove_cvref_t<T>;
    static constexpr std::size_t nn = N;
    static constexpr std::size_t ni = I;
    static constexpr std::size_t no = O;
    static constexpr std::size_t nr = R;

  public:
    using observer_t = MIMO_KalmanObserver<value_type, N, I, O>;

    struct parameter_t
    {
      observer_t::parameter_t const              observer;
      exmath::matrix_t<value_type, nr, nn> const C_int;
      exmath::matrix_t<value_type, ni, nn> const hT_kal;

      exmath::matrix_t<value_type, nr, nr> const          A_int;
      exmath::matrix_t<value_type, nr, 2 * ni + nr> const B_int;
      exmath::matrix_t<value_type, ni, nr> const          hT_int;
    };

    MIMO_KalmanIntegralController(parameter_t const& param)
        : m_param{ param }
    {
    }

    exmath::matrix_t<value_type, no, 1> calculate_y_kal() const { return this->m_kal.get_y(); }
    exmath::matrix_t<value_type, nr, 1> calculate_y_int() const { return this->m_param.C_int * this->m_kal.get_x(); }
    exmath::matrix_t<value_type, ni, 1> calculate_u_star() const
    {
      return exmath::matrix_t<value_type, ni, 1>{} - this->m_param.hT_kal * this->m_kal.get_x() - this->m_param.hT_int * this->m_x_int;
    }
    exmath::matrix_t<value_type, nn, 1> const& get_x() const { return this->m_kal.get_x(); }

    void correction_step(exmath::matrix_t<value_type, no, 1> const& y_meas) { return this->m_kal.correction_step(y_meas); }

    void update_observer(exmath::matrix_t<value_type, ni, 1> const& u) { return this->m_kal.update_step(u); }

    void update_integrator(exmath::matrix_t<value_type, ni, 1> const& u_star,
                           exmath::matrix_t<value_type, ni, 1> const& u_dash,
                           exmath::matrix_t<value_type, nr, 1> const& e)
    {
      auto const u = [](exmath::matrix_t<value_type, ni, 1> const& u_star, exmath::matrix_t<value_type, ni, 1> const& u_dash,
                        exmath::matrix_t<value_type, nr, 1> const& e) -> exmath::matrix_t<value_type, 2 * ni + nr, 1>
      {
        std::size_t                                  j = 0;
        exmath::matrix_t<value_type, 2 * ni + nr, 1> ret{};
        for (std::size_t i = 0; i < u_dash.number_of_rows; i++, j++)
        {
          ret(j, 0) = u_dash(i, 0);
        }
        for (std::size_t i = 0; i < u_star.number_of_rows; i++, j++)
        {
          ret(j, 0) = u_star(i, 0);
        }
        for (std::size_t i = 0; i < e.number_of_rows; i++, j++)
        {
          ret(j, 0) = e(i, 0);
        }
        return ret;
      }(u_star, u_dash, e);
      this->m_x_int = this->m_param.A_int * this->m_x_int + this->m_param.B_int * u;
    }

    void reset_observer() { this->m_kal.reset(); }

    void reset_integrator() { this->m_x_int = {}; }

  private:
    parameter_t const& m_param;

    observer_t                          m_kal{ this->m_param.observer };
    exmath::matrix_t<value_type, nr, 1> m_x_int = {};
  };

  template <typename T, std::size_t N> using SISO_KalmanIntegralController = MIMO_KalmanIntegralController<T, N, 1, 1, 1>;

}    // namespace controller::calculator

#endif
