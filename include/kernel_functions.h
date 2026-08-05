#ifndef KERNEL_FUNCTIONS
#define KERNEL_FUNCTIONS

#include <array>
#include <cmath>
#include <iostream>
#include <numbers>

template<size_t D>
struct cubic_spline
{
  static constexpr double a = [] {
    if constexpr (D==1) return 1.0;
    else if constexpr (D==2) return 15/(7*std::numbers::pi);
    else if constexpr (D==3) return 3/(2*std::numbers::pi);
    else
      static_assert(D<=3, "Dimension should be <=3");
  }();

  inline static double h_pow_D = 0;
  inline static double h_pow_DP1 = 0;

  cubic_spline(double h)
  {
    h_pow_D = std::pow(h,D);
    h_pow_DP1 = std::pow(h,D+1);
  }

  static double eval(double dist, double h)
  {
    double q = dist / h;
    if (q < 1)
    {
      return 1 / h_pow_D * a * (2. / 3. - q * q + 0.5 * q * q * q);
    }
    if (q < 2)
    {
      return 1 / h_pow_D * a * 1. / 6. * (2-q)*(2-q)*(2-q);
    }
    return 0;
  }

  static std::array<double, D> eval_derivatives(std::array<double, D> diff, double h)
  {
    std::array<double, D> res(diff);

    double dist = norm<D>(diff);
    if (dist == 0.)
      return {};

    double coef = a * 1. / h_pow_DP1 / dist;
    double q = dist / h;
    if (q < 1)
    {
      coef *= -2. * q + 3. / 2. * q * q;
    }
    else if (q < 2)
    {
      coef *= -0.5 * (2 - q) * (2 - q);
    }
    else
    {
      coef = 0;
    }

    for (size_t i = 0; i < D; ++i)
    {
      res[i] *= coef;
    }
    return res;
  }
};

#endif