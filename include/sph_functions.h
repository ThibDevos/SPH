#ifndef SPH_FUNCTIONS
#define SPH_FUNCTIONS

#include <constants.h>
#include <particle.h>

template<size_t D>
void compute_density_pressure(std::vector<Particle<D>> & particles)
{
  for(auto & p : particles)
  {
    p.density = 0;
    for(auto const  q : p.neighbours)
    {
      double dist = distance(p,*q);
      p.density += q->mass * cubic_spline<D>::eval(dist, H);
    }
    p.pressure = K * (std::pow(p.density/RHO_0, 7) -1);
  }
}

template<size_t D>
void compute_forces(std::vector<Particle<D>> & particles)
{
  for(auto & p : particles)
  {
    p.forces.fill(0.);
    double coef_i_pressure = -p.mass;
    double coef_i_viscosity = p.mass * NU;
    for(auto const  q : p.neighbours)
    {
      std::array<double, D> diff;
      for(int i=0; i<D; ++i){diff[i] = p.position[i] - q->position[i];}
      auto grad_kernel = cubic_spline<D>::eval_derivatives(diff, H);

      double coef_j_pressure =  q->mass * (p.pressure/(p.density*p.density) + q->pressure/(q->density*q->density));
      double coef_j_viscosity = 2. * q->mass / q->density;

      for(int i=0; i<D; ++i)
      {
        p.forces[i] += coef_i_pressure * coef_j_pressure * grad_kernel[i]; //pressure
        p.forces[i] += coef_i_viscosity * coef_j_viscosity * (p.velocity[i] - q->velocity[i]) * (inner<D>(diff, grad_kernel))/(inner<D>(diff, diff) + 0.01 * H * H);
      }
    }
    p.forces[1] -= p.mass * 9.81 ;
  }
}

template<size_t D>
void update(std::vector<Particle<D>> & particles)
{
  for(auto & p : particles)
  {
    for(int i=0; i<D; ++i)
    {
      p.velocity[i] += p.forces[i] * DT / p.mass;
      p.position[i] += p.velocity[i] * DT;
    }
  }
}

template <size_t D>
void check_boundaries(std::vector<Particle<D>> &particles, std::array<double, D> const &bot_left, std::array<double, D> const &top_right)
{
  for (auto &p : particles)
  {
    for (int i = 0; i < D; ++i)
    {
      if ( p.position[i] < bot_left[i])
      {
        p.position[i] = bot_left[i];
        p.velocity[i] *= -0.3; // rebond amorti, pas juste écrasé
      }
      else if ( p.position[i] > top_right[i])
      {
        p.position[i] = top_right[i];
        p.velocity[i] *= -0.3;
      }
    }
  }
}

#endif