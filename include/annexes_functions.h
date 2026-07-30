#ifndef ANNEXES_FUNCTIONS
#define ANNEXES_FUNCTIONS

#include <cmath>
#include <particle.h>

template<size_t D>
double distance(Particle<D> const & p, Particle<D> const & q)
{
  double res = 0.;
  for (std::size_t i = 0; i < D; ++i)
  {
    const double diff = p.position[i] - q.position[i];
    res += diff * diff;
  }
  return std::sqrt(res);
}

template<size_t D, typename T> 
inline double norm(T a)
{
  double res = 0;
  for(int i=0; i<D; ++i)
  {
    res += a[i]*a[i];
  }
  return std::sqrt(res);
}

template<size_t D, typename T> 
inline double inner(T a, T b)
{
  double res = 0;
  for(int i=0; i<D; ++i)
  {
    res += a[i]*b[i];
  }
  return res;
}


#endif