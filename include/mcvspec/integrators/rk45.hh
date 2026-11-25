#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <utility>

#include "array_operators.hh"
#include "integration.hh"

namespace mcvspec { namespace integration { namespace rk45 {

namespace tableau {

constexpr int order = 5;

constexpr std::array<std::array<double, 5>, 6> a = {{
    {           0,             0,            0,          0,             0},
    {         0.2,             0,            0,          0,             0},
    {       0.075,         0.225,            0,          0,             0},
    {     44./45.,      -56./15.,       32./9.,          0,             0},
    {19372./6561., -25360./2187., 64448./6561., -212./729.,             0},
    { 9017./3168.,     -355./33., 46732./5247.,   49./176., -5103./18656.}
}};

constexpr std::array<double,6> b = {35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.};

constexpr std::array<double,7> c = {0.,1./5.,3./10.,4./5.,8./9.,1.,1.};

constexpr std::array<double,7> e = {-71./57600., 0., 71./16695.,-71./1920., 17253./339200.,-22./525., 1./40.};

constexpr std::array<std::array<double,4>,7> p = {{
    {1.,-8048581381./2820520608., 8663915743./2820520608.,-12715105075./11282082432.},
    {0., 0., 0., 0.},
    {0., 131558114200./32700410799., -68118460800./10900136933., 87487479700./32700410799.},
    {0.,-1754552775./470086768., 14199869525./1410260304.,-10690763975./1880347072.},
    {0., 127303824393./49829197408.,-318862633887./49829197408., 701980252875./199316789632.},
    {0.,-282668133./205662961., 2019193451./616988883.,-1453857185./822651844.},
    {0., 40617522./29380423.,-110615467./29380423., 69997945./29380423.}
}};

}

template<DifferentialEquation System>
class RK45 {
public:
    using State = typename System::State;

    RK45(System system, const State& y0, const Options& options)
        : func(std::move(system)),
            y(y0),
            y_prev(y0),
            y_trial(y0),
            options(options)
    {
        assert(y0.size() > 0);
        for(auto& stage : k){
            stage = y0;
        }

        for(auto& coeff : q){
            coeff = y0;
        }
    }

    void reset(double t0, double tf, const State& y0);
    void reset(double t0, double tf, double dt0, const State& y0);
    void set_system(System eq);

    Status step();

    State interpolate(const double t) const;
    void interpolate(const double t, State& out) const;

    double coordinate() const noexcept { return t; };
    const State& state() const noexcept { return y; };
    const State& dydt() const noexcept { return k.back(); };
    void dydt(State& dydt) const noexcept { dydt = k.back(); };

private:
    double initial_step_size();
    void compute_step();
    void compute_interpolant() const;
    double error() const;

    System func;
    double t, t_prev, t_end, dt;
    State y{}, y_prev{};

    std::array<State,7> k;
    mutable std::array<State,tableau::p[0].size()> q;
    mutable bool interpolant_current = false;
    State y_trial{};

    Options options;
};

template<DifferentialEquation System>
double RK45<System>::initial_step_size() {
    const double span = std::abs(t_end - t);

    if (span == 0.0){
        dt=0;
        return 0;
    }

    const double direction = std::copysign(1.0, t_end-t);
    double d[3] = {0,0,0};
    for(size_t i=0; i<y.size(); ++i){
        const double s = options.atol + options.rtol*std::abs(y[i]);
        d[0] += y[i]*y[i]/(s*s);
        d[1] += k.back()[i]*k.back()[i]/(s*s);
    }
    d[0] = std::sqrt(d[0]/y.size());
    d[1] = std::sqrt(d[1]/y.size());

    double dt0 = (d[0] < 1e-5 || d[1] < 1e-5) ? span*1e-6 : 1e-2*d[0]/d[1];
    dt0 = std::min(dt0, span);

    linear_combo(1.0, y, direction*dt0, k.back(), y_trial);
    func(t + direction*dt0, y_trial, k[1]);

    for (size_t i = 0; i < y.size(); ++i) {
        const double s = options.atol + options.rtol*std::abs(y[i]);
        d[2] += (k[1][i] - k.back()[i])*(k[1][i] - k.back()[i])/(s*s);
    }

    d[2] = std::sqrt(d[2]/y.size())/dt0;

    double dt1 = std::pow(
        1e-2/std::max(d[1],d[2]),
        1.0/tableau::order
    );

    return direction*std::min({100.0*dt0, dt1, span});
}


template<DifferentialEquation System>
void RK45<System>::reset(double t0, double tf, const State& y0){
    assert(y0.size() == y.size());
    t=t0;
    t_end=tf;
    y=y0;
    func(t, y, k.back());
    interpolant_current=false;

    dt = initial_step_size();
}

template<DifferentialEquation System>
void RK45<System>::reset(double t0, double tf, double dt0, const State& y0){
    assert(y0.size() == y.size());
    t=t0;
    t_end=tf;
    y=y0;
    func(t, y, k.back());
    interpolant_current=false;

    dt = std::copysign(dt0, t_end-t);
}

template<DifferentialEquation System>
void RK45<System>::set_system(System system) {
    func = std::move(system);
    interpolant_current = false;
}


template<DifferentialEquation System>
void RK45<System>::compute_step(){
    // populates RK45::y_trial with an attempted step
    // AND RK45::k with RK stage derivatives
    for(std::size_t i=1; i<tableau::a.size(); ++i){
        y_trial = y;
        for(int j=0; j<i; ++j){
            linear_combo(y_trial, dt*tableau::a[i][j], k[j]);
        }
        func(t + tableau::c[i]*dt, y_trial, k[i]);
    }
    y_trial = y;
    for (std::size_t j=0; j < tableau::b.size(); ++j){
        linear_combo(y_trial, dt*tableau::b[j],k[j]);
    }
    func(t+dt, y_trial, k[tableau::b.size()]);
}

template<DifferentialEquation System>
double RK45<System>::error() const {
    double sqr_err = 0;
    for(size_t i=0; i<y.size(); ++i){
        const double tol = options.atol + options.rtol*std::abs(y_trial[i]);
        double err=0;
        for(size_t j=0; j<tableau::e.size(); ++j){
            err += tableau::e[j]*k[j][i];
        }
        err *= dt/tol;
        sqr_err += err*err;
    }
    return sqr_err/y.size();
}

template<DifferentialEquation System>
Status RK45<System>::step(){
    constexpr double safe_margin = 0.9;
    constexpr double min_step_factor = 0.2;
    constexpr double max_step_factor = 5.0;
    constexpr int max_attempts = 100;

    bool was_rejected = false;
    int attempts=0;

    k[0] = k.back();
    interpolant_current = false;

    while(attempts<max_attempts){
        ++attempts;
        if(!std::isfinite(dt)){ return Status::failure;}
        if(t+dt == t){ return Status::step_size_underflow;}

        if (std::abs(dt) > std::abs(t_end - t)){ dt = t_end - t; }
        compute_step();
        const double err = error();

        if(err <= 1){
            t_prev = t;
            y_prev = y;
            t += dt;
            y = y_trial;
            if(was_rejected){
                dt *= std::min(safe_margin*std::pow(err,-0.5/tableau::order), 1.0);
            }
            else if(err == 0){
                dt *= max_step_factor;
            }
            else{
                dt *= std::min(safe_margin*std::pow(err,-0.5/tableau::order), max_step_factor);
            }

            return Status::success;
        }

        was_rejected = true;
        dt *= std::max(safe_margin*std::pow(err,-0.5/tableau::order), min_step_factor);
    }

    return Status::max_steps_reached;
}
template<DifferentialEquation System>
void RK45<System>::compute_interpolant() const {
    q[0] = k[0];
    for(size_t i=1; i<q.size(); ++i){
        std::fill(q[i].begin(), q[i].end(), 0.0);
        for(size_t j=0; j<k.size(); ++j){
            linear_combo(q[i], tableau::p[j][i], k[j]);
        }
    }
    interpolant_current = true;
}

template<DifferentialEquation System>
void RK45<System>::interpolate(const double t_interp, State& out) const {
    if(!interpolant_current){ compute_interpolant(); }

    const double s = (t_interp - t_prev)/(t - t_prev);

    out = q[3];
    for(int i=2; i>=0; --i){
        linear_combo(1, q[i], s, out, out);
    }
    linear_combo(1, y_prev, (t - t_prev)*s, out, out);
}

template<DifferentialEquation System>
typename System::State RK45<System>::interpolate(const double t_interp) const {
    State out{};
    interpolate(t_interp, out);
    return out;
}

}}}
