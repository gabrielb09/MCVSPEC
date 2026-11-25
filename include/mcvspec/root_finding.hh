#pragma once

#include <concepts>

namespace mcvspec{

template<typename T>
concept ScalarFunction =
    requires(const T& func, double x) {
        { func(x) } -> std::same_as<double&>;
    };

namespace root_finding{

struct Bounds{
    double upper, lower;
};

struct Options{
    double tol;
};

}}
