#pragma once

#include <concepts>
#include <ranges>


namespace mcvspec {

template<typename T>
concept Container =
    std::ranges::random_access_range<T> &&
    std::ranges::sized_range<T> &&
    std::same_as<std::ranges::range_value_t<T>, double> &&
    std::default_initializable<T> &&
    std::copyable<T> &&
    requires(T& s, const T& cs, std::size_t i) {
        { s[i] } -> std::same_as<double&>;
        { cs[i] } -> std::same_as<const double&>;
    };

template<typename T>
concept DifferentialEquation =
    requires { typename T::State; } &&
    Container<typename T::State> &&
    requires(
        const T& system,
        double t,
        const typename T::State& y,
        typename T::State& dydx
    ) {
        { system(t, y, dydx) } -> std::same_as<void>;
    };

}
