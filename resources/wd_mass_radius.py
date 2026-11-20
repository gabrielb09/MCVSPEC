from collections.abc import Callable
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar


class OdeEvent:
    terminal : bool = False
    def __init__(self,
        func : Callable,
        term : bool = False
    ):
        self.func = func
        self.terminal = term
    def __call__(self,
        t : float,
        y : np.ndarray,
        args : tuple
    ) -> float:
        return self.func(t,y,args)


c = 29979245800.0 # cm/s
hbar = 6.62607015e-34/1.602176634e-19/(2*np.pi)/1e3 # keV s
m_e = 510.99895069 # keV
m_e_ergs = 9.1093837139e-28*c*c # g
alpha = 0.0072973525693
gell_mann_coeff = (1-np.log(2))/np.pi**2 # https://doi.org/10.1103/PhysRev.106.364
M_sun = 1.988409870698051e+33 # g
R_sun = 69570000000.0 # cm
amu = 1.66053906892e-24 # g
G = 6.67430e-8 # cm^3/g/s^2

r0 = 1e-32
_density_coeff = (amu/(3*np.pi**2))*(m_e/(hbar*c))**3 # g/cm3
pressure_coeff = (m_e_ergs/(3*np.pi**2))*(m_e/(hbar*c))**3 # erg/cm3


def pressure(
    x : float,
    Z : int
) -> float :
    x2 = x*x
    x3 = x2*x
    x4 = x3*x
    gamma2 = 1+x2
    gamma = np.sqrt(gamma2)
    eta = np.arcsinh(x)

    P_chandra = 0.125*(x*gamma*(2*x2 - 3) + 3*eta)
    tf_factor = alpha*np.cbrt(12*Z*Z/np.pi)
    P_coulomb = -tf_factor*x4*(0.1 + 6*tf_factor*x/(175*gamma))
    P_exchange = -0.75*(alpha/np.pi)*(1.5*(x2+eta*eta) - (x*eta/gamma)*(3 + x2) - x4/6)
    P_correlate = -(gell_mann_coeff/3)*alpha*alpha*x3
    return P_chandra + P_coulomb + P_exchange + P_correlate


def dP_dx(
    x : float,
    Z : int
) -> float :
    x2 = x*x
    x3 = x2*x
    x4 = x3*x
    gamma2 = 1+x2
    gamma = np.sqrt(gamma2)
    eta = np.arcsinh(x)

    dP_chandra = x4/gamma
    tf_factor = alpha*np.cbrt(12*Z*Z/np.pi)
    dP_coulomb = -0.4*tf_factor*x3*(1 + 3*tf_factor*x*(5 + 4*x2)/(35*gamma2*gamma))
    dP_exchange = (alpha/(2*np.pi))*(x3/gamma2)*(3*eta*x/gamma + x - 2)
    dP_correlate = -gell_mann_coeff*alpha*alpha*x2

    return dP_chandra + dP_coulomb + dP_exchange + dP_correlate

def equation_of_state(
    r : float,
    state : tuple[float, float],
    Z : int
) -> tuple[float, float]:

    x,m = state
    dx_dr = -m*(x**3)/(dP_dx(x, Z)*(r**2))
    dm_dr = 4*np.pi*(r**2)*(x**3)

    return dx_dr, dm_dr

def get_mass_and_radius(
    mu : float,
    Z : int,
    x0 : float,
    xf : float|None = None

) -> tuple[float, float]:
    if xf is None:
        xf = root_scalar(lambda x: pressure(x, Z),
            bracket=(1e-6, 0.2),
            x0=0.1,
            method="brentq",
            rtol=1e-10, xtol=1e-16
        ).root
    zero_pressure = OdeEvent(lambda r, state, Z : state[0]-xf, term=True)
    density_coeff =  mu*_density_coeff
    r_max = R_sun*density_coeff*np.sqrt(G/pressure_coeff)
    m0 = (4/3)*np.pi*(x0*r0)**3
    res = solve_ivp(
        equation_of_state,
        (r0, r_max),
        [x0,m0],
        args = (Z,),
        method = "DOP853",
        rtol = 1e-10,
        atol = 1e-16,
        events = [zero_pressure]
    )

    m_f = res.y_events[-1][-1,-1]
    r_f = res.t_events[-1][-1]

    mass = m_f*np.sqrt((pressure_coeff/G)**3)/density_coeff**2
    radius = r_f*np.sqrt(pressure_coeff/G)/density_coeff

    return mass/M_sun, radius/R_sun

def _x0_bounds(
    mass_radius : Callable,
    mass_range : tuple[float,float],
    max_r_e : float
) -> tuple[float, float]:
    """Finds the initial relativity parameter corresponding to two final mass
    bounds

    The upper bound is determined by geometrically searching outwards until an
    x0 is found that yields a larger mass than the target

    The Lower bound is determined by geometrically searching inwards until an
    x0 is found that yields a smaller mass or a minimum viable x0 is reached
    corresponding to r_e < (Z^(2/3))/10. For lower densities the degenerate
    fermi model is no longer valid.
    """

    x1 = np.cbrt(9*np.pi/4)*alpha/max_r_e
    if mass_radius(x1)[0]>mass_range[0]:
        x_min=x1
    else:
        x_min = root_scalar(lambda x: mass_radius(x)[0]-mass_range[0],
            method="brentq",
            bracket=(x1, 1e2),
            rtol=1e-10, xtol=1e-16
        ).root

    if x_min > 1:
        x1 = x_min
    else:
        x1 = 1
    m1,r1 = mass_radius(x1)
    while(True):
        x2 = x1*2
        m2,r2 = mass_radius(x2)

        if m1 < mass_range[1] < m2 :
            x_max = x2
            break

        dm_dr = (m2-m1)/(r2-r1)
        m_err = np.abs(m2 - (m2 - dm_dr*r2))
        if m_err/m2 < 1e-6:
            x_max = 2*x2
            break

        m1 = m2
        r1 = r2
        x1 = x2

    return x_min, x_max

def get_mass_radius_relationship(
    mu : float,
    Z : int,
    mass_range : tuple[float,float] = (0.2, 1.42),
    tol : float = 1e-6
) -> tuple[np.ndarray, np.ndarray]:

    xf = root_scalar(lambda x: pressure(x, Z),
        method="brentq",
        bracket=(1e-8, 0.2),
        x0=0.1,
        rtol=1e-10, xtol=1e-16
    ).root

    @lru_cache(maxsize=1024)
    def mass_radius(x0):
        return get_mass_and_radius(mu, Z, x0, xf)

    x = list(_x0_bounds(mass_radius, mass_range, 0.1/np.cbrt(Z*Z)))
    mass = [mass_radius(x0)[0] for x0 in x]
    radius = [mass_radius(x0)[1] for x0 in x]


    i = 0
    while i < len(x) - 1:
        x_mid = (x[i+1]+x[i])/2
        m,r = mass_radius(x_mid)
        dr_dm = (radius[i+1]-radius[i])/(mass[i+1]-mass[i])
        r_interp = radius[i] + dr_dm*(m-mass[i])
        err = np.abs(r_interp/r-1)

        if err < tol:
            i += 1
            continue

        x.insert(i+1,x_mid)
        mass.insert(i+1, m)
        radius.insert(i+1, r)

    return np.array(mass), np.array(radius)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="White Dwarf mass-radius Solver",
        description="Produces a list of masses and radii for a range of central densities"
    )
    parser.add_argument("-Z", type=float, required=True, help="charge of nucleus of WD material")
    parser.add_argument("-A", type=float, help="mass of nucleus of WD material (in amu), superceded by mu if input")
    parser.add_argument("-mu", type=float, help="A/Z ratio for WD material, if not specified A/Z is computed from other inputs")
    parser.add_argument("--header", action="store_true", help="generate header file instead of text")
    args = parser.parse_args()

    if args.mu is None:
        args.mu = args.A/args.Z

    mass, radius = get_mass_radius_relationship(args.mu, args.Z)

    if args.header:
        with open("mass_radius.hh", "w+") as out_file:
            out_file.write("#pragma once\n\n")
            out_file.write("#include <array>\n\n")
            out_file.write("namespace mcvspec{\n\n")
            out_file.write(f"constexpr std::array<double, {mass.size}> white_dwarf_mass{{")
            for i, m in enumerate(mass*M_sun):
                if i == 0:
                    out_file.write(f"\n    {m}")
                    continue
                if i%5 == 0:
                    out_file.write(f",\n    {m}")
                    continue
                out_file.write(f", {m}")
            out_file.write("\n};\n")
            out_file.write(f"constexpr std::array<double, {radius.size}> white_dwarf_radius{{")
            for i, r in enumerate(radius*R_sun):
                if i == 0:
                    out_file.write(f"\n    {r}")
                    continue
                if i%5 == 0:
                    out_file.write(f",\n    {r}")
                    continue
                out_file.write(f", {r}")
            out_file.write("\n};\n\n}")

    else:
        with open("mass_radius.txt", "w+") as out_file:
            out_file.write("Mass [Solar Mass], Radius [Solar Radii]\n")
            out_file.writelines(f"{m}, {r}\n" for m,r in zip(mass, radius))
