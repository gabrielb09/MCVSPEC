#include <array>
#include <iostream>

#include "integration.hh"
#include "rk45.hh"


class TestEq{
public:
    using State = std::array<double, 2>;
    explicit TestEq(const double k, const double m) :  spring_const(k), mass(m) {}

    void operator()(double t, const State& y, State& dydt) const {
        // y = x,v
        dydt[0] = y[1];
        dydt[1] = -(spring_const/mass)*y[0];
    }
private:
    double spring_const, mass;
};

int main(){
    using namespace mcvspec::integration;

    rk45::RK45 spring(
        TestEq(1,1),
        std::array<double, 2>{0,0},
        mcvspec::integration::Options{.rtol=1e-8, .atol=1e-10}
    );

    spring.reset(0, 10, {1,0});
    auto status = Status::success;
    while(spring.coordinate()<10 && status==Status::success){
        status = spring.step();
    }

    if (status != Status::success) {
        std::cerr << "Integration failed\n";
        return static_cast<int>(status);
    }

    const auto& y = spring.state();

    const double x_exact = std::cos(10.0);
    const double v_exact = -std::sin(10.0);


    std::cout << "Final time: " << spring.coordinate() << '\n';
    std::cout << "Position error: " << std::abs(y[0] - x_exact) << '\n';
    std::cout << "Velocity error: "<< std::abs(y[1] - v_exact) << '\n';

    return static_cast<int>(status);
}
