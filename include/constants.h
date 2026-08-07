#ifndef CONSTANTS
#define CONSTANTS

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sph
{
struct SimulationConstants
{
  std::size_t target_particle_count = 1000;
  std::size_t actual_particle_count = 1000;
  double dx = 0.00653;
  double H = 0.00816;
  double RHO_0 = 1000.0;
  double K = 40921.7142;
  double DT = 0.00005;
  double NU = 0.01;
  double MASS = 0.042632;
  double DUMPING_WALL = 0.0;
};

inline SimulationConstants make_simulation_constants(std::size_t particle_count,
                                                    double dam_width = 0.146,
                                                    double dam_height = 0.292,
                                                    double rho_0 = 1000.0,
                                                    double gravity = 9.81,
                                                    double neighbors_ratio = 1.25,
                                                    double monaghan_factor = 10.0,
                                                    double cfl_safety = 0.8,
                                                    double dt_margin = 1.1,
                                                    double nu = 0.001)
{
  SimulationConstants constants{};
  constants.target_particle_count = std::max<std::size_t>(1, particle_count);

  const double area = dam_width * dam_height;
  const double dx = std::sqrt(area / static_cast<double>(constants.target_particle_count));
  const std::size_t nx = std::max<std::size_t>(1, static_cast<std::size_t>(std::floor(dam_width / dx + 1.0)));
  const std::size_t ny = std::max<std::size_t>(1, static_cast<std::size_t>(std::floor(dam_height / dx + 1.0)));

  constants.dx = dx;
  constants.actual_particle_count = nx * ny;
  constants.H = neighbors_ratio * dx;
  constants.RHO_0 = rho_0;
  constants.MASS = rho_0 * dx * dx;

  const double c_sound = monaghan_factor * std::sqrt(gravity * dam_height);
  constants.K = (c_sound * c_sound * rho_0) / 7.0;

  const double dt_max = cfl_safety * constants.H / c_sound;
  constants.DT = dt_max / dt_margin;
  constants.NU = nu;
  constants.DUMPING_WALL = 0.9;
  return constants;
}

inline SimulationConstants& configuration()
{
  static SimulationConstants instance = make_simulation_constants(1000);
  return instance;
}

inline void configure(std::size_t particle_count)
{
  configuration() = make_simulation_constants(particle_count);
}

inline SimulationConstants const& constants()
{
  return configuration();
}
} // namespace sph

#endif