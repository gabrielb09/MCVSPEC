#pragma once

#include <concepts>
#include <cstddef>
#include <ranges>


namespace mcvspec {

template<typename T>
concept NumericalContainer =
    std::ranges::random_access_range<T> &&
    std::ranges::random_access_range<const T> &&
    std::ranges::sized_range<T> &&
    std::same_as<std::ranges::range_value_t<T>, double> &&
    std::default_initializable<T> &&
    std::copyable<T> &&
    !std::ranges::view<T> &&
    requires(T& arr, const T& const_arr, std::size_t i){
        { arr.size() } -> std::convertible_to<std::size_t>;
        { const_arr.size() } -> std::convertible_to<std::size_t>;
        { arr[i] } -> std::same_as<double&>;
        { const_arr[i] } -> std::same_as<const double&>;
    };

template<typename T>
concept DifferentialEquation =
    requires { typename T::State; } &&
    NumericalContainer<typename T::State> &&
    requires(
        const T& system,
        double t,
        const typename T::State& y,
        typename T::State& dydx
    ) {
        { system(t, y, dydx) } -> std::same_as<void>;
    };

namespace integration {

struct Options{
    double rtol = 1e-8;
    double atol = 1e-10;
    std::size_t max_steps = 100000;
};

enum class Status{
    success=0,
    failure=1,
    max_steps_reached=2,
    step_size_underflow=3
};

template<typename State>
struct Result{
    State final_state;
    Status status;
    std::size_t steps;
};

template<DifferentialEquation System>
auto integrate(
    const System& system,
    double t0,
    double t1,
    const typename System::State& y0,
    const Options& options
) -> Result<typename System::State>;


}}
