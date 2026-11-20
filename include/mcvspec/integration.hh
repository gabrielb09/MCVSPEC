#pragma once

#include <cstddef>

#include "ode.hh"


namespace mcvspec { namespace integration {

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
