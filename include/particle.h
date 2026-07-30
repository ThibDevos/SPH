#ifndef PARTICLE
#define PARTICLE

#include <iostream>

#include <array>


template<size_t D>
struct Particle
{
  std::array<double, D> position;
  std::array<double, D> velocity;
  std::array<double, D> forces;
  double pressure;
  double density;
  double mass = 0.4;

  // std::vector<Particle<D>> neighbour;
};

#endif