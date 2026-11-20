#include <iostream>

#include "Cataclysmic_Variable.hh"
#include "constants.hh"

int main(){
    using namespace mcvspec;
    White_Dwarf am_her = White_Dwarf{
        .mass = 0.62*m_sol,
        .radius = Mass_to_Radius(0.62*m_sol),
        .b_field = 13.6e6,
        .cos_inclination = 0.892,
        .inverse_mag_radius = 0,
        .corotation_radius = 0,
        .distance = 87.92*pc_to_cm
    };
    Accretion_Column col = Accretion_Column{
        .accretion_rate=1.27,
        .accretion_area=1e15,
        .metallicity=1,
        .shock_pressure_ratio=0.75,
        .sin_mag_colat=0.57
    };
    Tolerance tol = Tolerance{
        .absolute_error=1e-8,
        .relative_error=1e-6,
        .kT_grid_spacing=1,
        .altitude_grid_spacing=0.05
    };


    Cataclysmic_Variable mCV(am_her, col, tol);
    std::cout << "Finding shock height" << std::endl;
    mCV.Find_Shock_Position();
    std::cout << "Interpolating column" << std::endl;
    mCV.Build_Column_Profile();
    mCV.Print_Properties();

    return 0;
}
