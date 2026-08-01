#ifndef MORTON
#define MORTON

#include <algorithm>
#include <array>
#include <bitset>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "particle.h"

template <size_t D>
void radix_sort_indices(std::vector<uint64_t>& codes, std::vector<int>& indices, std::vector<Particle<D>>& particles)
{
  size_t n = codes.size();
  if (n == 0) return;

  std::vector<uint64_t> temp_codes(n);
  std::vector<int> temp_indices(n);
  std::vector<Particle<D>> temp_particles(n);

  for (int shift = 0; shift < 64; shift += 8)
    {
      size_t count[257] = {0};

      for (size_t i = 0; i < n; i++)
        {
          size_t bucket = (codes[i] >> shift) & 0xFF;
          count[bucket + 1]++;
        }

      for (size_t i = 0; i < 256; i++)
        {
          count[i + 1] += count[i];
        }

      for (size_t i = 0; i < n; i++)
        {
          size_t bucket = (codes[i] >> shift) & 0xFF;
          size_t pos = count[bucket]++;
          temp_codes[pos] = codes[i];
          temp_indices[pos] = indices[i];
          temp_particles[pos] = std::move(particles[i]);
        }

      codes.swap(temp_codes);
      indices.swap(temp_indices);
      particles.swap(temp_particles);
    }
}

template<size_t D>
class Morton
{
  public:
    Morton(std::vector<Particle<D>> & particles, double h)
    {
      int nb_part = particles.size();
      indices.resize(nb_part);
      Morton_code.resize(nb_part);
      std::iota(indices.begin(), indices.end(), 0);
      std::vector<std::array<uint32_t,D>> int_coords(nb_part);

      //transform the floating points coordiantes of the vertices to Morton code
      float_coord_to_int(particles, int_coords, h);
      int_coord_to_Morton(int_coords);

      radix_sort_indices(Morton_code, indices, particles);

      build_leaves(particles);
      compute_neighbours(particles);
    }
  
    //XXX also compute l
    void float_coord_to_int(std::vector<Particle<D>> & particles, std::vector<std::array<uint32_t,D>>& int_coords, double h )
    {
      constexpr uint64_t max_coord = ((uint64_t(1)<<L) - 1);
      double min[D];
      double max[D];
      for(int i=0; i<D; ++i)
      {
        min[i] = std::numeric_limits<double>::max();
        max[i] = std::numeric_limits<double>::lowest();
      }
      for(auto const & p : particles)
      {
        for(int d = 0; d<D; ++d)
        {
          double pos = p.position[d];
          if(pos<min[d]) min[d] = pos;
          if(pos>max[d]) max[d] = pos;
        }
      }
      // Max length
      double extent[D];
      for(int d=0; d<D; ++d)
      {
        extent[d] = max[d] - min[d];
      }

      double Max = std::max(extent[0], extent[1]);
      if constexpr (D==3) Max = std::max(Max, extent[2]);
      if (Max == 0.0)
        Max = 1.0;

      {
        int i = 0;
        for (auto const &p : particles)
        {
          for (int d = 0; d < D; ++d)
          {
            int_coords[i][d] = static_cast<uint32_t>(std::clamp((p.position[d] - min[d])/Max, 0., 1.) * max_coord);
          }
          ++i;
        }
      }
      //use max and min to determine l
      l = (int)std::floor(std::log(Max / (2*h)) / std::log(2));
      std::cout << "l = " << l << std::endl;
    }

    void int_coord_to_Morton(std::vector<std::array<uint32_t, D>> &int_coords)
    {
      std::cout << "int_coordinates_to_Morton" << std::endl;
      int i = 0;
      for (auto &c : int_coords)
      {
        Morton_code[i] = encodeMorton(c);
        ++i;
      }
    }
    /*
    This methods builds the leaves.
    A leaf as a Morton code written using l<=L bits.
    A vertex, is in a leaf if its l first bits (from the left) are the same (example 1011 is in leaf 10)
    The first leaf is created from the first vertex (by taking its l first bits). Then,
    if the next vertex as the same first l bits of the leaf, we store it. Otherwise, we browsed all vertices which
    are in the leaf. We determine then the next leaf, using the first l bits of the first vertex that is not in the previous leaf.

    The l first bits are taken using   >> (D*(L-l))
    */
    void build_leaves(std::vector<Particle<D>> & particles)
    {
      for (int d = 0; d < D; ++d)
      {
        grid_limits[d] = 0;
      }

      uint64_t leaf_code = Morton_code[0] >> (D * (L - l)); // Morton code of the first leaf
      leaf current_leaf;
      current_leaf.Morton_code = leaf_code;
      current_leaf.vertex_indices.push_back(0);
      particles[0].leaf = leaf_code;

      for (size_t i = 1; i < Morton_code.size(); ++i)
      {
        uint64_t next_vertex_code = Morton_code[i] >> (D * (L - l));
        if (next_vertex_code == leaf_code) // Morton_code[i] is still in current_leaf
        {
          current_leaf.vertex_indices.push_back(i);
          particles[i].leaf = leaf_code;
        }
        else // Morton_code[i] is not in current_leaf. We create a new leaf from it
        {
          uint32_t x[D];
          DecodeMorton(current_leaf.Morton_code, x);
          for (int d = 0; d < D; ++d)
            if (x[d] > grid_limits[d])
              grid_limits[d] = x[d];

          leaves.push_back(current_leaf);

          current_leaf = leaf();
          current_leaf.Morton_code = next_vertex_code;
          current_leaf.vertex_indices.push_back(i);
          leaf_code = next_vertex_code;
          particles[i].leaf = leaf_code;
        }
      }
      // we add the last leaf
      uint32_t x[D];
      DecodeMorton(current_leaf.Morton_code, x);
      for (int d = 0; d < D; ++d)
        if (x[d] > grid_limits[d])
          grid_limits[d] = x[d];
      leaves.push_back(current_leaf);
      std::cout<<leaves.size()<<" leaves\n";
    }

    void compute_neighbours(std::vector<Particle<D>> & particles)
    {
      int constexpr pow3 = std::pow(3,D);
      uint32_t x[D];
      int nb_n = 0;
      for (auto it = leaves.begin(); it != leaves.end(); ++it)
      {
        auto &current_leaf = *it;
        DecodeMorton(current_leaf.Morton_code, x); // need to go back to interger coordinates
        // set if the leaf is near a boundary or not
        int boundary = 0;
        for(int d=0; d<D; ++d)
        {
          boundary |= (x[d] == 0);
          boundary <<= 1;
        }
        if constexpr (D==2) boundary <<= 1;
        for(int d=0; d<D; ++d)
        {
          boundary |= (x[d] == grid_limits[d]);
          boundary <<= 1;
        }
        if constexpr (D==2) boundary <<= 1;


        current_leaf.boundary = boundary;



        std::vector<Particle<D>*> neighbouring_parts;
        neighbouring_parts.reserve(it->vertex_indices.size() * pow3);
        uint32_t nx[D];
        for(int c=0; c<pow3; ++c) //3^D is the number of leaves in the cluster
        {
          int offset[D];
          int tmp = c;

          bool is_current = true;
          for(int d=0; d<D; ++d)
          {
            int digit = tmp%3;
            offset[d] = digit - 1; //-1, 0, 1
            if(offset[d]!=0) is_current = false;
            tmp/=3;
          }
          if(is_current) continue;

          bool out_of_limits = false;
          for(int d=0; d<D; ++d)
          {
            if(x[d] + offset[d] <0 || x[d] + offset[d] > grid_limits[d]) 
            {
              out_of_limits = true;
              break;
            }
            nx[d] = x[d] + offset[d];
          }
          if(out_of_limits) continue;
          leaf search_leaf;
          uint64_t n_code = encodeMorton(nx); // recode Morton
          search_leaf.Morton_code = n_code;
          auto found_it = std::lower_bound(leaves.begin(), leaves.end(), search_leaf, [](const leaf &a, const leaf &b)
                                           { return a.Morton_code < b.Morton_code; }); // we use the fact that the leaves are ordered

          if (found_it != leaves.end() && found_it->Morton_code == n_code)
          {
            current_leaf.neighbours.push_back(&leaves[found_it - leaves.begin()]);
            for(int i=0; i<found_it->vertex_indices.size(); ++i)
            {
              neighbouring_parts.push_back(&particles[found_it->vertex_indices[i]]);
            }
          }
        }

        for(int i=0; i<current_leaf.vertex_indices.size(); ++i)
        {
          particles[current_leaf.vertex_indices[i]].neighbours.clear();
          particles[current_leaf.vertex_indices[i]].neighbours.insert(particles[current_leaf.vertex_indices[i]].neighbours.begin(), neighbouring_parts.begin(), neighbouring_parts.end());
          for(int j=0; j<current_leaf.vertex_indices.size(); ++j)
          {
            particles[current_leaf.vertex_indices[i]].neighbours.push_back(&particles[current_leaf.vertex_indices[j]]);
          }
          
        }
      }
    }

// Spread the bits of a coordinate so that they can be interleaved with the
// other coordinates. The 2D and 3D variants use different masks.
uint64_t expand_bits(uint32_t x)
    {
      if constexpr (D == 2)
      {
        uint64_t y = x;
        y = (y | (y << 16)) & 0x0000FFFF0000FFFFULL;
        y = (y | (y << 8))  & 0x00FF00FF00FF00FFULL;
        y = (y | (y << 4))  & 0x0F0F0F0F0F0F0F0FULL;
        y = (y | (y << 2))  & 0x3333333333333333ULL;
        y = (y | (y << 1))  & 0x5555555555555555ULL;
        return y;
      }
      else if constexpr (D == 3)
      {
        uint64_t y = x & 0x1fffff;
        y = (y | y << 32) & 0x1f00000000ffff;
        y = (y | y << 16) & 0x1f0000ff0000ff;
        y = (y | y << 8) & 0x100f00f00f00f00f;
        y = (y | y << 4) & 0x10c30c30c30c30c3;
        y = (y | y << 2) & 0x1249249249249249;
        return y;
      }
      else
      {
        throw std::runtime_error("Morton encoding is only supported for 2D or 3D");
      }
    }

    // Inverse of expand_bits.
uint32_t compact_bits(uint64_t x)
{
    if constexpr(D == 2)
    {
        x &= 0x5555555555555555ULL;

        x = (x | (x >> 1))  & 0x3333333333333333ULL;
        x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0FULL;
        x = (x | (x >> 4))  & 0x00FF00FF00FF00FFULL;
        x = (x | (x >> 8))  & 0x0000FFFF0000FFFFULL;
        x = (x | (x >> 16)) & 0x00000000FFFFFFFFULL;

        return static_cast<uint32_t>(x);
    }

    else if constexpr(D == 3)
    {
        x &= 0x1249249249249249ULL;

        x = (x | (x >> 2)) & 0x10c30c30c30c30c3ULL;
        x = (x | (x >> 4)) & 0x100f00f00f00f00fULL;
        x = (x | (x >> 8)) & 0x1f0000ff0000ffULL;
        x = (x | (x >> 16)) & 0x1f00000000ffffULL;
        x = (x | (x >> 32)) & 0x1fffffULL;

        return static_cast<uint32_t>(x);
    }

    else
    {
        throw std::runtime_error("Morton decoding is only supported for 2D or 3D");
    }
}

    // from integer coordinates to Morton code
    template<typename T> //T is a container of D unint32_t
    uint64_t encodeMorton(T x)
    {
      uint64_t a = 0;
      for(int d=0; d<D; ++d)
        {
          a |= expand_bits(x[d]) << d;
        }
      return a;
    }

    // from a Morton code to integer coordinates
    template<typename T> //T is a container of D unint32_t
    void DecodeMorton(uint64_t code, T x)
    {
      for(int d=0; d<D; ++d)
        x[d] = compact_bits(code >> d);
    }

  private:
    static constexpr uint32_t L = []{
      if constexpr(D==2) return 32;
      else if constexpr (D==3) return 21;
      else std::runtime_error("Wrong dimension");
    }(); //each coordinate is represented by 32 bits in 2D (2*32=64 bits) or 21 bits in 3D (3*21 = 63)
    int l=2;  //the cube will be devided in 2^l*2^l(*2^l) cubes (called leaves) in which the particles are


    std::vector<uint64_t> Morton_code; // Morton code is a way to express in 1D, 2 or 3D coordinates (integer)
    std::vector<int> indices; //Morton_code will be sorted. indices helps keeping a relation of indices between Morton_code and particles

    uint32_t grid_limits[D]; //upper bounds of the grid (lower bounds are 0)

    struct leaf
    {
      uint64_t Morton_code; //coordinates of the leaf
      std::vector<int> vertex_indices; //indices of the vertices in it (indices correspond to the indices of the mesh vertices)
      std::vector<leaf*> neighbours; //pointers to the neighbour leaves (26 in 3d)
      int boundary = 0; // boundary \in [0,63], such that, each bit correspond to a boundary and is 1 if the leaf is close to it.
      //  The bits correspond, from right to left to x_-, y_-, z_-, x_+, y_+, z_+
    };

    std::vector<leaf> leaves;
};

#endif