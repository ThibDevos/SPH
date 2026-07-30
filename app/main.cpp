
#include <particle.h>
#include <initialisation.h>
#include <annexes_functions.h>
#include <kernel_functions.h>
#include <sph_functions.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>


void write_particles_vtu(const std::vector<Particle<2>>& particles,
                        const std::string& filename)
{
  std::ofstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  file << "  <PolyData>\n";
  file << "    <Piece NumberOfPoints=\"" << particles.size() << "\" NumberOfVerts=\""
       << particles.size() << "\">\n";

  file << "      <Points>\n";
  file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.position[0] << ' ' << particle.position[1]
         << " 0.0\n";
  }
  file << "        </DataArray>\n";
  file << "      </Points>\n";

  file << "      <PointData>\n";
  file << "        <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.pressure << '\n';
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.density << '\n';
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.mass << '\n';
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.velocity[0] << ' ' << particle.velocity[1]
         << " 0.0\n";
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.forces[0] << ' ' << particle.forces[1]
         << " 0.0\n";
  }
  file << "        </DataArray>\n";
  file << "      </PointData>\n";

  file << "      <Verts>\n";
  file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < particles.size(); ++i) {
    if (i > 0) {
      file << ' ';
    }
    file << i;
  }
  file << "\n        </DataArray>\n";

  file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 1; i <= particles.size(); ++i) {
    if (i > 1) {
      file << ' ';
    }
    file << i;
  }
  file << "\n        </DataArray>\n";
  file << "      </Verts>\n";

  file << "    </Piece>\n";
  file << "  </PolyData>\n";
  file << "</VTKFile>\n";
}

int main()
{
  std::vector<Particle<2>> particles(1000);

  std::array<double, 2> dam_bot_left{0.0, 0.0};
  std::array<double, 2> dam_top_right{0.5, 0.8};
  double dx = 0.02;

  std::array<double, 2> dom_bot_left{0.0, 0.0};
  std::array<double, 2> dom_top_right{2.0, 1.0};

  init_grid_block(dam_bot_left, dam_top_right, dx, particles);
  // search_neighbour(particles);

  double T = 10.;
  double t = 0.0;
  int count = 0;
  int count_files = 0;
  int save_freq = 100;
  while (t < T)
  {
    std::cout << t << "/" << T << std::endl;
    compute_density_pressure(particles);
    compute_forces(particles);
    update(particles);
    check_boundaries(particles, dom_bot_left, dom_top_right);
    if (count % save_freq == 0)
    {

      const std::filesystem::path output_path =
          std::filesystem::path(__FILE__).parent_path().parent_path() / "bin" / (std::string("particles_") + std::to_string(count_files) + ".vtp");
      count_files++;
      write_particles_vtu(particles, output_path.string());
    }
    ++count;
    t += DT;
  }
}