#pragma once
#ifndef CONTROLLER_INTERFACE_HPP_INCLUDED
#define CONTROLLER_INTERFACE_HPP_INCLUDED

#include <concepts>
#include <cstdint>
#include <type_traits>

namespace controller
{
  struct Status_Interface
  {
    enum class State : uint32_t
    {
      disabled = 0,
      enabled  = 1,
      fault,
    };

    enum class Fault : uint32_t
    {
      no_error = static_cast<uint32_t>(State::enabled),
      generic  = static_cast<uint32_t>(State::fault),
      invalid_measure_value,
      plant_missmatch,
    };

    struct status_t
    {
    public:
      constexpr status_t(State state)
          : m_value{ static_cast<uint32_t>(state) }
      {
      }

      constexpr status_t(Fault fault)
          : m_value{ static_cast<uint32_t>(fault) }
      {
      }

      constexpr status_t(uint32_t value)
          : m_value{ value }
      {
      }

      constexpr State as_state() const
      {
        if (this->m_value > static_cast<uint32_t>(State::enabled))
          return State::fault;
        return State{ this->m_value };
      }

      constexpr Fault as_fault() const
      {
        if (this->m_value <= static_cast<uint32_t>(Fault::no_error))
          return Fault::no_error;
        return Fault{ this->m_value };
      }
      constexpr uint32_t as_value() const { return this->m_value; }
      constexpr bool     has_fault() const { return this->as_fault() != Fault::no_error; }
      constexpr          operator State() const { return as_state(); }
      constexpr          operator Fault() const { return as_fault(); }
      constexpr          operator uint32_t() const { return as_value(); }

    private:
      uint32_t m_value;
    };

    virtual ~Status_Interface() = default;

    [[nodiscard]] virtual status_t get_status() const = 0;
  };

  template <typename T> struct Actuator_Interface: public Status_Interface
  {
    using target_type = std::remove_cv_t<T>;

    virtual ~Actuator_Interface() = default;

    struct actuator_output_t
    {
      target_type ranged_value;
      target_type target_value;
    };

    [[nodiscard]] virtual actuator_output_t set_target_value(target_type const& wish_value) = 0;
    [[nodiscard]] virtual target_type       get_target_value() const                        = 0;
  };

  template <typename TT, typename TV> struct Controller_Interface: public Actuator_Interface<TT>
  {
    using value_type = std::remove_cv_t<TV>;

    virtual ~Controller_Interface() = default;

    [[nodiscard]] virtual value_type get_value() const = 0;
  };
}    // namespace controller

#endif
