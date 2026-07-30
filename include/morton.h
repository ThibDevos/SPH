#ifndef MORTON
#define MORTON

#include <iostream>
#include <vector>

void radix_sort_indices(std::vector<uint64_t>& codes, std::vector<int>& indices)
{
  size_t n = codes.size();
  if (n == 0) return;

  std::vector<uint64_t> temp_codes(n);
  std::vector<int> temp_indices(n);

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
        }

      codes.swap(temp_codes);
      indices.swap(temp_indices);
    }
}

template<size_t D, typename T>
class Morton
{
  public:
    Morton(std::vector<T> & particles, double h)
    {
      int nb_part = particles.size();
      indices.resize(nb_part);
      Morton_code.resize(nb_part);
      std::iota(indices.begin(), indices.end(), 0);
      std::vector<std::array<uint32_t,D>> int_coords(nb_part);

      //transform the floating points coordiantes of the vertices to Morton code
      float_coord_to_int(particles, int_coords, h);
      int_coord_to_Morton(int_coords);

      radix_sort_indices(Morton_code, indices);

      build_leaves
    }
  
    //XXX also compute l
    void float_coord_to_int(std::vector<T> & particles, std::vector<std::array<uint32_t,D>> int_coords, double h )
    {
      constexpr uint32_t max_coord = ((1u<<L) - 1);
      double min[D];
      double max[D];
      for(int i=0; i<D; ++i)
      {
        min[i] = std::numeric_limits<double>::max();
        max[i] = std::numeric_limits<double>::min();
      }
      for(auto const & p : particles)
      {
        for(int d = 0; d<D; ++d)
        {
          double pos = p.position[i];
          if(pos<min[d]) min[d] = pos;
          if(pos>max[d]) max[d] = pos;
        }
      }
      // Max length
      double extent[3];
      extent[0] = max[0] - min[0];
      extent[1] = max[1] - min[1];
      extent[2] = max[2] - min[2];

      double Max = std::max(extent[0], std::max(extent[1], extent[2]));
      if (Max == 0.0)
        Max = 1.0;

      {
        int i = 0;
        for (auto const &p : particles)
        {
          for (int d = 0; d < D; ++d)
          {
            int_coord[i][d] = static_cast<uint32_t>(std::clamp((p.position[d] - min[d])/Max, 0., 1.) * max_coord);
          }
          ++i;
        }
      }
      //use max and min to determine l
      l = (int)std::floor(std::log(Max / h) / std::log(2));
      std::cout << "l = " << l << std::endl;
    }

    void int_coordinates_to_Morton(std::vector<std::array<uint32_t, 3>> &int_coords)
    {
      std::cout << "int_coordinates_to_Morton" << std::endl;
      int i = 0;
      for (auto &c : int_coords)
      {
        Morton_code[i] = encodeMorton(c[0], c[1], c[2]);
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
    void build_leaves()
    {
      uint64_t leaf_code = Morton_code[0] >> (D*(L-l)); //Morton code of the first leaf
    }



    // x is written using 21 bits (all 9 bits forming the 32bits int are 0).
    // This function returns the same number but with two 0 between each bit
    // ex : x =1101 -> y = 1001000001
    //  see https://stackoverflow.com/questions/18529057/produce-interleaving-bit-patterns-morton-keys-for-3d-coordinates-for-32-bit
    //  see https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
    uint64_t Morton_Linear_Octree::expand_bits(uint32_t x)
    {
      uint64_t y = x & 0x1fffff;
      y = (y | y << 32) & 0x1f00000000ffff;  // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
      y = (y | y << 16) & 0x1f0000ff0000ff;  // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
      y = (y | y << 8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
      y = (y | y << 4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
      y = (y | y << 2) & 0x1249249249249249;
      return y;
    }

    // inverse of expand_bits.
    uint32_t Morton_Linear_Octree::compact_bits(uint64_t x)
    {
      x &= 0x1249249249249249;
      x = (x ^ (x >> 2)) & 0x10c30c30c30c30c3;
      x = (x ^ (x >> 4)) & 0x100f00f00f00f00f;
      x = (x ^ (x >> 8)) & 0x1f0000ff0000ff;
      x = (x ^ (x >> 16)) & 0x1f00000000ffff;
      x = (x ^ (x >> 32)) & 0x1fffff;
      return static_cast<uint32_t>(x);
    }

    // from integer coordinates to Morton code
    uint64_t Morton_Linear_Octree::encodeMorton(uint32_t x, uint32_t y, uint32_t z)
    {
      uint64_t a = 0;
      a |= expand_bits(x) | expand_bits(y) << 1 | expand_bits(z) << 2;
      return a;
    }

    // from a Morton code to integer coordinates
    void Morton_Linear_Octree::DecodeMorton(uint64_t code, uint32_t &x, uint32_t &y, uint32_t &z)
    {
      x = compact_bits(code >> 0);
      y = compact_bits(code >> 1);
      z = compact_bits(code >> 2);
    }

  private:
    static constexpr uint32_t L = []{
      if constexpr(D==2) return 32;
      else if constexpr (D==3) return 21;
      else runtime_error("Wrond dimension");
    }(); //each coordinate is represented by 32 bits in 2D (2*32=64 bits) or 21 bits in 3D (3*21 = 63)
    int l=5;  //the cube will be devided in 2^l*2^l(*2^l) cubes (called leaves) in which the particles are


    std::vector<uint64_t> Morton_code; // Morton code is a way to express in 1D, 2 or 3D coordinates (integer)
    std::vector<int> indices; //Morton_code will be sorted. indices helps keeping a relation of indices between Morton_code and particles
};

#endif