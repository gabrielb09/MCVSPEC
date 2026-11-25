#pragma once

#include <cmath>

#include "root_finding.hh"

namespace mcvspec{ namespace root_finding{
namespace itp {

// note, these are defined here for consistency with the ITP parameterization
// however, my implementation of the ITP method assumes these values and is
// simplified accordingly. As such, changing these values has no impact on the
// algorithm as written.

constexpr double k1 = 0.2;
constexpr double k2 = 2;
constexpr double n0 = 1;

template<ScalarFunction Func>
inline double find_root(Func& func, Bounds& bounds, const Options& opts){
    Bounds fvals{
        .upper=func(bounds.upper),
        .lower=func(bounds.lower)
    };
    return find_root(func, bounds, fvals, opts);
}

template<ScalarFunction Func>
inline double find_root(Func& func, Bounds& bounds, Bounds& fvals, const Options& opts){
    double r0 = bounds.upper-bounds.lower;

    while(bounds.upper-bounds.lower > 2*opts.tol){
        const double bisection = 0.5*(bounds.upper+bounds.lower);
        const double regula_falsi = (fvals.upper*bounds.lower - fvals.lower*bounds.upper)/(fvals.upper - fvals.lower);
        const double diff = bisection-regula_falsi;
        const double dir = std::copysign(1,diff);

        const double delta = 0.2*(bounds.upper-bounds.lower);
        const double r = r0 - 0.5*(bounds.upper-bounds.lower);

        double new_bound = 0;
        if(std::abs(diff) <= delta){ new_bound = bisection; }
        else if(std::abs(diff) - delta <= r){ new_bound = regula_falsi + dir*delta; }
        else{ new_bound = bisection - dir*r; }

        r0 *= 0.5;

        const double new_val = func(new_bound);
        if(new_val == 0.0){
            return new_bound;
        }
        else if(std::signbit(new_val) == std::signbit(fvals.upper)){
            bounds.upper = new_bound;
            fvals.upper = new_val;
        }
        else{
            bounds.lower = new_bound;
            fvals.lower = new_val;
        }
    }
    return (bounds.upper+bounds.lower)/2;
}

}}}
