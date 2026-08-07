#ifndef INITIALISATION
#define INITIALISATION

#include <particle.h>

#include <cmath>
#include <vector>

template<size_t D>
void equaly_distributed(std::array<double, D>& bot_left, std::array<double, D>& top_right, std::vector<Particle<D>>& particles);

template<size_t D>
void init_grid_block(std::array<double, D> const& bot_left,
                      std::array<double, D> const& top_right,
                      double dx,
                      std::vector<Particle<D>>& particles)
{
  std::array<std::size_t, D> counts;
  std::size_t total = 1;
  for (std::size_t d = 0; d < D; ++d)
  {
    double span = top_right[d] - bot_left[d];
    counts[d] = std::max<std::size_t>(1, static_cast<std::size_t>(std::floor(span / dx)) + 1);
    total *= counts[d];
  }

  for (std::size_t j = 0; j < total; ++j)
  {
    std::array<std::size_t, D> index{};
    std::size_t stride = 1;
    for (std::size_t d = D; d-- > 0;)
    {
      index[d] = (j / stride) % counts[d];
      stride *= counts[d];
    }
    for (std::size_t d = 0; d < D; ++d)
      particles[j].position[d] = bot_left[d] + dx * static_cast<double>(index[d]);
  }
}

template <size_t D>
void init_grid_block(std::array<double, D> const& bot_left,
                      std::array<double, D> const& top_right,
                      double dx,
                      Particles<D>& particles)
{
  std::array<std::size_t, D> counts;
  std::size_t total = 1;
  for (std::size_t d = 0; d < D; ++d)
  {
    double span = top_right[d] - bot_left[d];
    counts[d] = std::max<std::size_t>(1, static_cast<std::size_t>(std::floor(span / dx)) + 1);
    total *= counts[d];
  }

  for (std::size_t j = 0; j < total; ++j)
  {
    std::array<std::size_t, D> index{};
    std::size_t stride = 1;
    for (std::size_t d = D; d-- > 0;)
    {
      index[d] = (j / stride) % counts[d];
      stride *= counts[d];
    }
    for (std::size_t d = 0; d < D; ++d)
      particles.position(j,d) = bot_left[d] + dx * static_cast<double>(index[d]);
  }
}
#endif