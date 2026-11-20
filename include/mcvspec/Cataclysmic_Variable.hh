#pragma once

#include "dipole.hh"
#include "gaunt.hh"
#include "mass_radius.hh"
#include "rk45.hh"
#include <vector>

namespace mcvspec{


struct White_Dwarf{
    double mass;
    double radius;
    double b_field;
    double cos_inclination;
    double inverse_mag_radius;
    double corotation_radius;
    double distance;
};

struct Accretion_Column{
    double accretion_rate;
    double accretion_area;
    double metallicity;
    double shock_pressure_ratio;
    double sin_mag_colat;
};

struct Tolerance{
    double absolute_error;
    double relative_error;
    double kT_grid_spacing;
    double altitude_grid_spacing;
};

double Mass_to_Radius(double);
double Luminosity_to_Accretion_Rate(double, White_Dwarf);

struct Accretion_Coeffs{
    double cyclotron_const, bremss_const, exchange_const, coulomb_log_coeff;
    double avg_atomic_charge, mass_to_number_density;
    double density_conv, pressure_conv;
};

class Accretion_Eq{
public:
    using State = std::array<double, 4>;

    explicit Accretion_Eq(const Accretion_Coeffs& c, const Dipole& g) : coeffs(c), geometry(g) {}

    void operator()(double s, const State& y, State& dyds) const {
        const double& w = y[0];
        const double& x = y[1];
        const double& v = y[2];
        const double& pe = y[3];
        geometry.update_coordinates(w, r, proj_r_w, convergance, scale_factors);
        const double& h = scale_factors[0];
        const double area = scale_factors[1]*scale_factors[2];
        const double mdot = 1./area;
        const double p = mdot*(x-v);

        const double ne_cgs = coeffs.mass_to_number_density*coeffs.density_conv*mdot/v;
        const double kT_cgs = coeffs.pressure_conv*pe/ne_cgs;
        const double gff = gaunt::gaunt_factor(kT_cgs);
        const double coulomb_log = 0.5*log(coeffs.coulomb_log_coeff*kT_cgs*kT_cgs/ne_cgs);

        const double cyc = (coeffs.cyclotron_const/gff)*pe*pe*v*v*v*h*pow(v*v/h, 0.425);
        const double rad = coeffs.bremss_const*gff*sqrt(pe/(h*h*h*v*v*v))*(1+cyc);
        const double exch = (coeffs.exchange_const*coulomb_log/(coeffs.bremss_const*gff))*(coeffs.avg_atomic_charge*(p-pe)/pe - 1.0)/(pe*v*h*(1+cyc)); // ratio of exchange to radiation
        const double grav = 0.5*proj_r_w/(r*r)/rad;
        const double conv = v*convergance/(h*rad);

        dyds[0] = -1.5*p*v/(s*h*rad);
        dyds[1] = -1.5*(p/s)*(p*conv/mdot - grav);
        dyds[2] = 1.5*(p*v/s)*(2.0/mdot + 5*p*conv/mdot -3*grav)/(5*x - 8*v);
        dyds[3] = (p/s)*(1.0 - exch + 5*pe*(1.5*grav - 1.0/mdot - v*conv)/(5*x - 8*v));
    }

private:
    Accretion_Coeffs coeffs;
    Dipole geometry;
    mutable double r, proj_r_w, convergance, scale_factors[3];
};

class Cataclysmic_Variable{
    protected:
        static constexpr size_t n_dim = 4;
        static constexpr size_t n_grid_vars = 3;
        // input properties
        const White_Dwarf white_dwarf;
        const Accretion_Column accretion_column;
        const Tolerance error_control;
        const Dipole geometry;
        // unit conversion (cgs values of MCVSPEC nd unit system)
        const double length_conv, vel_conv, accretion_rate_conv;
        const double time_conv, mass_conv, volume_conv, energy_conv, density_conv, pressure_conv;

        // derived column properties
        std::vector<double> abundances; // fractional abundance of elements in accretion column
        double avg_ion_mass, avg_atomic_charge, mass_to_number_density, exchange_const, bremss_const, cyclotron_const;
        // thermal profile
        std::vector<double> position, altitude, volume_element, velocity, density, total_pressure, electron_pressure, electron_density, electron_temperature, ion_temperature;
        // shock boundary condition
        std::array<double, n_dim> shock_boundary;
        double shock_entropy;
        // if solution found
        bool valid_solution = true;


    public:
        Cataclysmic_Variable(White_Dwarf, Accretion_Column, Tolerance);
        virtual ~Cataclysmic_Variable() = default;

        void Find_Shock_Position();
        void Build_Column_Profile();

        void Print_Properties() const;

    protected:
        void Set_Abundances(){
            abundances.resize(n_elements);
            abundances = {1.00e+00, 9.77e-02, 3.63e-04, 1.12e-04, 8.51e-04, 1.23e-04,
                          3.80e-05, 2.95e-06, 3.55e-05, 1.62e-05, 3.63e-06, 2.29e-06,
                          4.68e-05, 1.78e-06}; // taken from Anders & Grevesse (1989) DOI: 10.1016/0016-7037(89)90286-X
            double total= abundances[0]+abundances[1];
            for(size_t i=2; i<abundances.size(); i++){
                abundances[i] *= accretion_column.metallicity;
                total += abundances[i];
            }
            std::transform(abundances.begin(),abundances.end(),abundances.begin(),[total](double x) {return x/total;});
            Set_Cooling_Constants();
        }
        void Set_Cooling_Constants();
        int Bracket_Shock_Position(double&,double&,double&,double&) const;
        template <typename func>
        void Build_Grid(func,const std::array<double, n_grid_vars>&, std::vector<std::array<double,n_dim>>&) const;
        void Compute_Shock_Bound(double, std::array<double,n_dim>&, double&) const;
        double Landing_Altitude(double) const;
        mutable integration::rk45::RK45<Accretion_Eq> integrator;
};

inline double Luminosity_to_Accretion_Rate(double luminosity, White_Dwarf wd){
    double accretion_rate = luminosity/(grav_const*wd.mass/wd.radius);
    return accretion_rate;
}

inline double Mass_to_Radius(double mass){
    int left_ind = 0;
    int i = white_dwarf_mass.size()/2;
    int right_ind = white_dwarf_mass.size()-1;
    while(right_ind-left_ind > 1){
        i = left_ind + (right_ind-left_ind)/2;
        if(mass>white_dwarf_mass[i]){
            left_ind = i;
        }
        else{
            right_ind = i;
        }
    }
    double delta_r = white_dwarf_radius[right_ind]-white_dwarf_radius[left_ind];
    double delta_m = white_dwarf_mass[right_ind]-white_dwarf_mass[left_ind];
    double radius = white_dwarf_radius[left_ind] + (delta_r/delta_m)*(mass-white_dwarf_mass[left_ind]);
    return radius;
}

}
