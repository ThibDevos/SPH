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
    for(auto const&  q : p.neighbours)
    {
      double dist = distance(p,*q);
      p.density += q->mass * cubic_spline<D>::eval(dist, sph::constants().H);
    }
    p.pressure = sph::constants().K * (std::pow(p.density / sph::constants().RHO_0, 7) - 1);
    if (p.pressure < 0.0) p.pressure = 0.0;
  }
}

template<size_t D>
void compute_density_pressure(Particles<D> & particles)
{
  for(int i=0; i<particles.nb_particles; ++i)
  {
    particles.density(i) = 0;
    for(auto j : particles.neighbours[i])
    {
      double dist = 0.;
      for(int d=0; d<D; ++d)
      {
        double dx = particles.position(i, d) - particles.position(j,d);
        dist += dx * dx;
      }
      dist = std::sqrt(dist);

      particles.density(i) += particles.mass * cubic_spline<D>::eval(dist, sph::constants().H);
    }
    particles.pressure(i) = sph::constants().K * (std::pow(particles.density(i) / sph::constants().RHO_0, 7) - 1);
    if (particles.pressure(i) < 0.0) particles.pressure(i) = 0.0;
  }
}

template<size_t D>
void compute_forces(std::vector<Particle<D>> & particles)
{
  for(auto & p : particles)
  {
    p.forces.fill(0.);
    double coef_i_pressure = -p.mass;
    double coef_i_viscosity = p.mass * sph::constants().NU;
    for(auto const &  q : p.neighbours)
    {
      std::array<double, D> diff;
      for(int i=0; i<D; ++i){diff[i] = p.position[i] - q->position[i];}
      auto grad_kernel = cubic_spline<D>::eval_derivatives(diff, sph::constants().H);

      double coef_j_pressure = q->mass * (p.pressure / (p.density * p.density) + q->pressure / (q->density * q->density));
      double coef_j_viscosity = 2. * q->mass / q->density;

      for(int i=0; i<D; ++i)
      {
        p.forces[i] += coef_i_pressure * coef_j_pressure * grad_kernel[i]; //pressure
        p.forces[i] += coef_i_viscosity * coef_j_viscosity * (p.velocity[i] - q->velocity[i]) * (inner<D>(diff, grad_kernel)) / (inner<D>(diff, diff) + 0.01 * sph::constants().H * sph::constants().H);
      }
    }
    p.forces[1] -= p.mass * 9.81 ;
  }
}

template <size_t D>
void compute_forces(Particles<D>& particles)
{
  for(int i=0; i<particles.nb_particles; ++i)
  {
    for(int d=0; d<D; ++d)
    {
      particles.force(i,d) = 0;
    }
    double coef_i_pressure = -particles.mass;
    double coef_i_viscosity = particles.mass * sph::constants().NU;

    for(auto j : particles.neighbours[i])
    {
      std::array<double, D> diff;
      for(int d=0; d<D; ++d){diff[d] = particles.position(i,d) - particles.position(j,d);}
      auto grad_kernel = cubic_spline<D>::eval_derivatives(diff, sph::constants().H);

      double coef_j_pressure = particles.mass * (particles.pressure(i) / (particles.density(i) * particles.density(i)) + particles.pressure(j) / (particles.density(j) * particles.density(j)));
      double coef_j_viscosity  = 2. * particles.mass / particles.density(j);

      for(int d=0; d<D; ++d)
      {
        particles.force(i,d) += coef_i_pressure * coef_j_pressure * grad_kernel[d]; //pressure
        particles.force(i,d) += coef_i_viscosity * coef_j_viscosity * (particles.velocity(i,d) - particles.velocity(j,d)) * (inner<D>(diff, grad_kernel)) / (inner<D>(diff, diff) + 0.01 * sph::constants().H * sph::constants().H);
      }
    }
    particles.force(i,1) -= particles.mass * 9.81;
  }
}

template<size_t D>
void update(std::vector<Particle<D>> & particles)
{
  for(auto & p : particles)
  {
    for(int i=0; i<D; ++i)
    {
      p.velocity[i] += p.forces[i] * sph::constants().DT / p.mass;
      p.position[i] += p.velocity[i] * sph::constants().DT;
    }
  }
}

template<size_t D>
void update(Particles<D> & particles)
{
  for(int i=0; i<particles.nb_particles; ++i)
  {
    for(int d=0; d<D; ++d)
    {
      particles.velocity(i,d) += particles.force(i,d) * sph::constants().DT / particles.mass;
      particles.position(i,d) += particles.velocity(i,d) * sph::constants().DT;
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
        p.velocity[i] *= -sph::constants().DUMPING_WALL; // rebond amorti, pas juste écrasé
      }
      else if ( p.position[i] > top_right[i])
      {
        p.position[i] = top_right[i];
        p.velocity[i] *= -sph::constants().DUMPING_WALL;
      }
    }
  }
}

template<size_t D>
void check_boundaries(Particles<D> & particles, std::array<double, D> const &bot_left, std::array<double, D> const &top_right)
{
  for (int i=0; i<particles.nb_particles; ++i)
  {
    for (int d = 0; d < D; ++d)
    {
      if ( particles.position(i,d) < bot_left[d])
      {
        particles.position(i,d) = bot_left[d];
        particles.velocity(i,d) *= -sph::constants().DUMPING_WALL; // rebond amorti, pas juste écrasé
      }
      else if ( particles.position(i,d) > top_right[d])
      {
        particles.position(i,d) = top_right[d];
        particles.velocity(i,d) *= -sph::constants().DUMPING_WALL;
      }
    }
  }
}


#endif