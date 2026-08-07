#ifndef PARTICLE
#define PARTICLE

#include <constants.h>

#include <iostream>

#include <array>
#include <vector>

#define SCALAR_ACCESSOR(name, member) \
  inline auto& name(size_t i) {return member[i];} \
  inline const auto& name(size_t i) const {return member[i];}

#define VECTOR_ACCESSOR(name, member) \
  inline auto& name(size_t i, size_t d) {return member[i * D + d];} \
  inline const auto& name(size_t i, size_t d) const {return member[i * D + d];}

template <size_t D>
struct Particles
{
  std::vector<double> positions;
  std::vector<double> velocities;
  std::vector<double> forces;
  std::vector<double> pressures;
  std::vector<double> densities;
  std::vector<double> leaves;
  double mass = sph::constants().MASS;

  std::vector<std::vector<int>> neighbours;

  size_t nb_particles;

  Particles(size_t N)
  {
    positions.resize(N*D);
    velocities.resize(N*D);
    forces.resize(N*D);
    pressures.resize(N);
    densities.resize(N);
    leaves.resize(N);
    neighbours.resize(N);

    nb_particles = N;
  }
  VECTOR_ACCESSOR(position, positions);
  VECTOR_ACCESSOR(velocity, velocities);
  VECTOR_ACCESSOR(force, forces);
  SCALAR_ACCESSOR(pressure, pressures);
  SCALAR_ACCESSOR(density, densities);
  SCALAR_ACCESSOR(leaf, leaves);
  SCALAR_ACCESSOR(neighbour, neighbours);

  void swap(Particles & other)
  {
    std::swap(positions, other.positions);
    std::swap(velocities, other.velocities);
    std::swap(forces, other.forces);
    std::swap(pressures, other.pressures);
    std::swap(densities, other.densities);
    std::swap(leaves, other.leaves);
    std::swap(neighbours, other.neighbours);
    std::swap(mass, other.mass);
    std::swap(nb_particles, other.nb_particles);
  }
};

template<size_t D>
struct Particle
{
  std::array<double, D> position;
  std::array<double, D> velocity;
  std::array<double, D> forces;
  double pressure;
  double density;
  double mass = sph::constants().MASS;

  std::vector<Particle<D>*> neighbours;
  int boundary = 0; // boundary \in [0,63], such that, each bit correspond to a boundary and is 1 if the leaf is close to it.
                    //  The bits correspond, from right to left to x_-, y_-, z_-, x_+, y_+, z_+

  int leaf=-1;
};

#endif