#ifndef PARTICLE
#define PARTICLE

#include <constants.h>

#include <iostream>

#include <array>
#include <vector>


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