#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>


namespace mcvspec {

Cataclysmic_Variable::Cataclysmic_Variable(White_Dwarf wd, Accretion_Column col, Tolerance tol):
    white_dwarf(wd),
    accretion_column(col),
    error_control(tol),
    geometry(accretion_column.sin_mag_colat*accretion_column.sin_mag_colat),
    length_conv(white_dwarf.radius),
    vel_conv(sqrt(2*grav_const*white_dwarf.mass/white_dwarf.radius)),
    accretion_rate_conv(geometry.a_0*col.accretion_rate),
    time_conv(length_conv/vel_conv),
    mass_conv(accretion_rate_conv*length_conv*length_conv*time_conv),
    volume_conv(length_conv*length_conv*length_conv),
    energy_conv(mass_conv*vel_conv*vel_conv),
    density_conv(mass_conv/volume_conv),
    pressure_conv(energy_conv/volume_conv),
    integrator(
        Accretion_Eq(
            Accretion_Coeffs(1,1,1,1,1,1,density_conv,pressure_conv),
            geometry
        ),
        std::array<double, 4>{},
        integration::Options{.rtol=error_control.relative_error, .atol=error_control.absolute_error}
    )
{Set_Abundances();}

void Cataclysmic_Variable::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
    // constants related to column composition
    avg_ion_mass = 0;
    avg_atomic_charge = 0;
    double avg_charge_squared = 0;
    double avg_charge_sqr_over_mass = 0;

    for(size_t i=0; i<abundances.size(); i++){
        avg_ion_mass += abundances[i]*atomic_mass[i]*amu_to_g;
        avg_atomic_charge += abundances[i]*atomic_charge[i];
        avg_charge_squared += abundances[i]*atomic_charge[i]*atomic_charge[i];
        avg_charge_sqr_over_mass += abundances[i]*atomic_charge[i]*atomic_charge[i]/(atomic_mass[i]*amu_to_g);
    }

    mass_to_number_density = avg_atomic_charge/(avg_ion_mass + avg_atomic_charge*m_e);
    const double j_e = mass_to_number_density*accretion_column.accretion_rate*geometry.a_0;
    const double b_sqr = white_dwarf.b_field*white_dwarf.b_field;
    bremss_const = bremss_coeff*(avg_charge_squared/avg_atomic_charge)*sqrt(j_e*j_e*j_e);
    bremss_const *= sqrt(time_conv*time_conv*time_conv*mass_conv)*length_conv/energy_conv;
    cyclotron_const = cyclotron_coeff*(avg_atomic_charge/avg_charge_squared);
    cyclotron_const *= b_sqr*pow(b_sqr/(geometry.a_0*accretion_column.accretion_area), 0.425)/(j_e*j_e);
    cyclotron_const *= pow(geometry.a_0*geometry.a_0/j_e, 1.85);
    cyclotron_const *= pressure_conv*pressure_conv*pow(vel_conv,3.85);
    exchange_const = exchange_coeff*avg_charge_sqr_over_mass*sqrt(j_e*j_e*j_e*j_e*j_e);
    exchange_const *= length_conv*length_conv*time_conv*time_conv*time_conv*sqrt(time_conv/energy_conv)/energy_conv;

    Accretion_Coeffs coeffs{
        .cyclotron_const=cyclotron_const,
        .bremss_const=bremss_const,
        .exchange_const=exchange_const,
        .coulomb_log_coeff=coulomb_log_coeff,
        .avg_atomic_charge=avg_atomic_charge,
        .mass_to_number_density=mass_to_number_density,
        .density_conv=density_conv,
        .pressure_conv=pressure_conv,
    };

    integrator.set_system(Accretion_Eq(coeffs, geometry));
}

void Cataclysmic_Variable::Compute_Shock_Bound(double w_s, std::array<double, 4>& bound, double& s_s) const{
    double r_s, proj_r_w, convergance, scale_factors[3];
    geometry.update_coordinates(w_s, r_s, proj_r_w, convergance, scale_factors);
    double mdot = 1./(scale_factors[1]*scale_factors[2]);
    double vff = sqrt(1./r_s - length_conv*white_dwarf.inverse_mag_radius);
    bound[0] = w_s;
    bound[1] = vff;
    bound[2] = vff*0.25;
    bound[3] = (accretion_column.shock_pressure_ratio/(1+accretion_column.shock_pressure_ratio))*mdot*0.75*vff;
    s_s = mdot*0.75*vff*pow(0.25*vff/mdot, 5./3.);
}

double Cataclysmic_Variable::Landing_Altitude(double w_s) const {
    // return signed distance from WD surface in w
    double s_s;
    std::array<double, 4> y;
    Compute_Shock_Bound(w_s, y, s_s);
    integration::Status status = integration::Status::success;
    integrator.reset(s_s, 0, y);
    while(integrator.coordinate()/s_s > 1e-2 && status==integration::Status::success){
         status=integrator.step();
    }
    if(status != integration::Status::success){ return -1; }

    double error = 1;
    std::array<double, 4> slope{};
    integrator.dydt(slope);

    double s_prev = integrator.coordinate();
    double dw_prev = slope[0];

    while(error > error_control.absolute_error){
        integrator.step();
        integrator.dydt(slope);
        error = 0.5*std::abs((slope[0]-dw_prev)/(integrator.coordinate()-s_prev))*integrator.coordinate()*integrator.coordinate(); //difference between linear and quadratic extrapolation on w
        s_prev = integrator.coordinate();
        dw_prev = slope[0];
    }
    integrator.step();
    integrator.dydt(slope);
    double landing = integrator.state()[0] - slope[0]*integrator.coordinate();
    return landing - geometry.w_0;
}

// Bracket the shock height by minimizing the landing coordinate
// (which is maximizing the landing altitude since w ~ -r close to the surface and at the pole)
// exit if
// 1. a point is found with w_l < w_0, in which case we have bracketed our solution
// 2. a minimum is found with w_l > w_0 in which case no solution exists
int Cataclysmic_Variable::Bracket_Shock_Position(double& upper_bound, double& lower_bound, double& upper_landing, double& lower_landing) const {
    upper_bound = geometry.w_0;
    upper_landing = Landing_Altitude(upper_bound);
    if(upper_landing==-2){
        return 1;
    }

    double r=1.;
    double dr = 0.01;
    auto r_to_w = [this](double r){ return sqrt(1-geometry.u*r)/(r*r); };
    double samples[3] = {Landing_Altitude(r_to_w(r-dr)),
                        upper_landing,
                        Landing_Altitude(r_to_w(r+dr))};

    double dwl_drs[2] = {(samples[1]-samples[0])/dr, (samples[2]-samples[1])/dr};
    double step = 0;
    while(dwl_drs[0]*dwl_drs[1] > 0){ // while minima not bounded
        if(samples[2]<0){
            lower_bound = r_to_w(r+dr);
            lower_landing = samples[2];
            return 0;
        }
        upper_bound = r_to_w(r+dr);
        upper_landing = samples[2];
        step = 0.5*dr*(samples[0]-samples[2])/(samples[0]-2*samples[1]+samples[2]);
        r += std::min(1.,step);
        samples[0] = Landing_Altitude(r_to_w(r-dr));
        samples[1] = Landing_Altitude(r_to_w(r));
        samples[2] = Landing_Altitude(r_to_w(r+dr));
        dwl_drs[0] = (samples[1]-samples[0])/dr;
        dwl_drs[1] = (samples[2]-samples[1])/dr;
    }

    if(samples[1] > 0){ // if minima > 0
        return 1;
    }

    lower_bound = r_to_w(r);
    lower_landing = samples[1];
    return 0;
}

void Cataclysmic_Variable::Find_Shock_Position(){
    double upper_bound, lower_bound;
    double upper_landing, lower_landing;
    int success = Bracket_Shock_Position(upper_bound, lower_bound, upper_landing, lower_landing);
    // if no minimum skip
    if(success != 0){
        valid_solution = false;
        return;
    }
    const double k1 = 0.2/(upper_bound-lower_bound);
    const double n0 = 1;
    double nmax = log2((upper_bound-lower_bound)/(2*error_control.absolute_error)) + n0;
    int i=0;
    double new_bound, new_altitude, midpoint, regula_falsi, truncation, projection, dir;
    while(upper_bound-lower_bound > error_control.absolute_error){
        midpoint = (upper_bound+lower_bound)/2;
        regula_falsi = (upper_landing*lower_bound - lower_landing*upper_bound)/(lower_landing-upper_landing);
        dir = (0. < (midpoint-regula_falsi)) - ((midpoint-regula_falsi) < 0.);
        truncation = k1*(upper_bound-lower_bound)*(upper_bound-lower_bound); // k2 = 2

        if(truncation <= std::abs(midpoint-regula_falsi)){
            new_bound = regula_falsi + dir*truncation;
        }
        else{
            new_bound = midpoint;
        }

        projection = error_control.absolute_error*(pow(2,nmax-i)) - (upper_bound-lower_bound)/2;
        if(std::abs(new_bound-midpoint) > projection){
            new_bound = midpoint - dir*projection;
        }

        new_altitude = Landing_Altitude(new_bound);
        if(new_altitude>0){
            upper_bound = new_bound;
            upper_landing = new_altitude;
        }
        else if(new_altitude<0){
            lower_bound = new_bound;
            lower_landing = new_altitude;
        }
        else{
            upper_bound = new_bound;
            lower_bound = new_bound;
        }
        i++;
    }
    std::array<double, 4> bound{};
    double s_s;
    Compute_Shock_Bound(0.5*(upper_bound+lower_bound), bound, s_s);
    shock_boundary = bound;
    shock_entropy = s_s;
}

template <typename func>
void Cataclysmic_Variable::Build_Grid(func grid_func, const std::array<double,n_grid_vars>& grid_spacing, std::vector<std::array<double, 4>>& grid) const {
    if(!valid_solution){
        return;
    }
    constexpr int n_segments=16;
    // integration variables
    double t = shock_entropy;
    std::array<double, 4> y = shock_boundary;
    integrator.reset(t, 0, y);

    // grid variables
    std::array<double,n_grid_vars> grid_vars{};
    grid_func(t,y,grid_vars);
    grid.reserve(size_t(std::abs(grid_vars[0]/grid_spacing[0]))+1);
    grid.push_back(y);

    // segment variables
    double t_left, t_right, t_cross;
    std::array<double, 4> y_left, y_right;
    std::array<double,n_grid_vars> grid_left, grid_right;
    y_right = y;
    grid_right = grid_vars;

    auto crossing_found = [&grid_spacing, &grid_vars, &grid_left, &grid_right, &t_cross](const double t_left, const double t_right){
        double delta_t, distance_right, distance_left;
        double min_dt = t_right-t_left;
        bool found=false;
        for(size_t i=0; i<n_grid_vars; ++i){
            distance_right = std::abs((grid_right[i]-grid_vars[i])/grid_spacing[i]);
            distance_left = std::abs((grid_left[i]-grid_vars[i])/grid_spacing[i]);
            delta_t = (t_right-t_left)*(1-distance_left)/(distance_right-distance_left);
            if(distance_right>1 && std::abs(delta_t)<std::abs(min_dt)){
                found=true;
                min_dt = delta_t;
            }
        }
        if(!found){
            return false;
        }
        t_cross = t_left + min_dt;
        return true;
    };

    while(grid_vars[0] > grid_spacing[0] && integrator.coordinate() > error_control.absolute_error){
        double t0 = integrator.coordinate();
        const auto status = integrator.step();
        if (status != integration::Status::success) {
            throw std::runtime_error("RK45 integration failed");
        }
        double dt = (integrator.coordinate()-t0)/double(n_segments);
        for(size_t seg=0; seg<n_segments; seg++){
            t_left = t0 + dt*seg;
            t_right = t0 + dt*(seg+1);
            y_left = y_right;
            integrator.interpolate(t_right, y_right);
            grid_left = grid_right;
            grid_func(t_right, y_right, grid_right);
            while(crossing_found(t_left, t_right)){
                t_left = t_cross;
                integrator.interpolate(t_left, y_left);
                grid_func(t_left, y_left, grid_left);
                grid.push_back(y_left);
                grid_vars = grid_left;
            }
        }
    }
}

void Cataclysmic_Variable::Build_Column_Profile(){
    if(!valid_solution){
        return;
    }
    const double dkTe = (error_control.kT_grid_spacing/erg_to_kev)*mass_to_number_density/(vel_conv*vel_conv);
    const double dkTi = dkTe/avg_atomic_charge;
    double r, proj, conv, scale_factors[3];
    geometry.update_coordinates(shock_boundary[0], r, proj, conv, scale_factors);
    const std::array<double,n_grid_vars> grid_spacing = {dkTi, dkTe, error_control.altitude_grid_spacing*(r-1)};

    auto kT_e = [](const std::array<double,4>& y, double scale_factors[3]){
        return y[2]*y[3]*scale_factors[0]*scale_factors[2];
    };
    auto kT_i = [](const std::array<double,4>& y, double scale_factors[3]){
        return y[2]*(y[1]-y[2]-y[3]*scale_factors[0]*scale_factors[2]);
    };
    auto grid_func = [this, kT_e, kT_i, &r, &proj, &conv, &scale_factors](const double& t, const std::array<double,n_dim>& y, std::array<double,n_grid_vars>& vars){
        geometry.update_coordinates(y[0], r, proj, conv, scale_factors);
        vars[0] = kT_i(y, scale_factors);
        vars[1] = kT_e(y, scale_factors);
        vars[2] = r;
    };

    std::vector<std::array<double,n_dim>> grid;
    Build_Grid(grid_func, grid_spacing, grid);

    int n_points = grid.size();
    position.resize(n_points);
    volume_element.resize(n_points);
    velocity.resize(n_points);
    altitude.resize(n_points);
    total_pressure.resize(n_points);
    electron_pressure.resize(n_points);
    electron_density.resize(n_points);
    density.resize(n_points);
    electron_temperature.resize(n_points);
    ion_temperature.resize(n_points);

    double mdot;

    for(size_t i=0; i<n_points; i++){
        geometry.update_coordinates(grid[i][0], r, proj, conv, scale_factors);
        position[i] = grid[i][0];
        volume_element[i] = (accretion_column.accretion_area*length_conv/geometry.a_0)*scale_factors[0]*scale_factors[1]*scale_factors[2];
        altitude[i] = length_conv*(r-1);
        velocity[i] = vel_conv*grid[i][2];
        mdot = 1./(scale_factors[1]*scale_factors[2]);
        total_pressure[i] = pressure_conv*mdot*(grid[i][1]-grid[i][2]);
        electron_pressure[i] = pressure_conv*grid[i][3];
        density[i] = density_conv*mdot/grid[i][2];
        electron_density[i] = mass_to_number_density*density[i];
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*(total_pressure[i]-electron_pressure[i])/(electron_density[i]/avg_atomic_charge);
    }
}

void Cataclysmic_Variable::Print_Properties() const{
    using std::cout;
    using std::endl;
    if(altitude.size() < 1){
        return;
    }
    cout << "===================================================" << endl;
    cout << "                   mCV Properties                  " << endl;
    cout << "===================================================" << endl;
    cout << " mass:               " << white_dwarf.mass/m_sol << " M_solar" << endl;
    cout << " radius:             " << white_dwarf.radius/r_sol << " R_solar" << endl;
    cout << " B_field:            " << white_dwarf.b_field/1e6 << " MG" << endl;
    if(white_dwarf.inverse_mag_radius != 0){
        cout << " R_m/R:              " << (1./white_dwarf.inverse_mag_radius)/white_dwarf.radius << endl;
    }
    cout << " accretion rate:     " << accretion_column.accretion_rate*accretion_column.accretion_area << " g/s" << endl;
    cout << " accretion rate:     " << density[0]*velocity[0] << " --> " <<  density[density.size()-1]*velocity[velocity.size()-1] << " g/cm2/s" << endl;
    cout << " shock height:       " << altitude[0]/white_dwarf.radius << " (h/R_wd)" << endl;
    cout << " shock height:       " << altitude[0] << " cm" << endl;
    cout << " shock temperature:  " << ion_temperature[0] << " keV" << endl;
    cout << " density:            " <<  electron_density[0] << " --> " << electron_density[electron_density.size()-1] <<  " e-/cm3" << endl;
}

}
