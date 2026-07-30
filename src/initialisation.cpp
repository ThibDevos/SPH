#include <initialisation.h>

#include <array>
#include <cmath>
#include <cstddef>

template<size_t D>
void equaly_distributed(std::array<double, D>& bot_left,
                        std::array<double, D>& top_right,
                        std::vector<Particle<D>>& particles)
{
  if (particles.empty()) {
    return;
  }

  std::array<std::size_t, D> counts{};
  counts.fill(1);

  std::size_t total_points = 1;
  while (total_points < particles.size()) {
    std::size_t dim = 0;
    for (std::size_t d = 1; d < D; ++d) {
      if (counts[d] < counts[dim]) {
        dim = d;
      }
    }

    ++counts[dim];
    total_points = 1;
    for (std::size_t d = 0; d < D; ++d) {
      total_points *= counts[d];
    }
  }

  for (std::size_t j = 0; j < particles.size(); ++j) {
    std::array<std::size_t, D> index{};
    std::size_t stride = 1;
    for (std::size_t d = D; d-- > 0;) {
      index[d] = (j / stride) % counts[d];
      stride *= counts[d];
    }

    for (std::size_t d = 0; d < D; ++d) {
      if (counts[d] > 1) {
        const double span = top_right[d] - bot_left[d];
        const double step = span / static_cast<double>(counts[d] - 1);
        particles[j].position[d] = bot_left[d] + step * static_cast<double>(index[d]);
      } else {
        particles[j].position[d] = bot_left[d];
      }
    }
  }
}

template void equaly_distributed<2>(std::array<double, 2>& bot_left,
                                    std::array<double, 2>& top_right,
                                    std::vector<Particle<2>>& particles);