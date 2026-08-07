#ifndef WRITE_H
#define WRITE_H

#include <particle.h>
#include <fstream>
#include <iomanip>

template <size_t D>
void write_particles_vtu(const std::vector<Particle<D>>& particles,
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
    file << "          " << std::setprecision(17) << particle.position[0] << ' ' << particle.position[1];
    if constexpr (D==3)
    {
      file << particle.position[2]<<"\n";
    }
    else
    {
      file << " 0.0\n";
    }
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

  file << "        <DataArray type=\"Float64\" Name=\"leaf\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.leaf << '\n';
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.velocity[0] << ' ' << particle.velocity[1];
    if constexpr (D==3)
    {
      file << particle.velocity[2]<<"\n";
    }
    else
    {
      file << " 0.0\n";
    }
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.forces[0] << ' ' << particle.forces[1];
    if constexpr (D==3)
    {
      file << particle.forces[2]<<"\n";
    }
    else
    {
      file << " 0.0\n";
    }
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

template <size_t D>
void write_particles_vtu(const Particles<D>& particles,
                        const std::string& filename)
{
  std::ofstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  file << "  <PolyData>\n";
  file << "    <Piece NumberOfPoints=\"" << particles.nb_particles << "\" NumberOfVerts=\""
       << particles.nb_particles << "\">\n";

  file << "      <Points>\n";
  file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int i=0; i<particles.nb_particles; ++i)
  {
    file << "          " << std::setprecision(17) << particles.positions[i*D + 0] << ' ' << particles.positions[i*D + 1];
    if constexpr (D==3)
    {
      file << ' ' << particles.positions[i*D + 1] <<"\n";
    }
    else
    {
      file <<" 0.0\n";
    }
  }
  file << "        </DataArray>\n";
  file << "      </Points>\n";

  file << "      <PointData>\n";
  file << "        <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">\n";
  for(int i=0; i<particles.nb_particles; ++i)
  {
    file << "          " << std::setprecision(17) << particles.pressures[i]<<"\n";
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n";
  for(int i=0; i<particles.nb_particles; ++i)
  {
    file << "          " << std::setprecision(17) << particles.densities[i]<<"\n";
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\">\n";
  for(int i=0; i<particles.nb_particles; ++i)
  {
    file << "          " << std::setprecision(17) << particles.mass<<"\n";
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"leaf\" format=\"ascii\">\n";
  for(int i=0; i<particles.nb_particles; ++i)
  {
    file << "          " << std::setprecision(17) << particles.leaves[i]<<"\n";
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int i=0; i<particles.nb_particles; ++i)
  {
    file << "          " << std::setprecision(17) << particles.velocities[i*D + 0] << ' ' << particles.velocities[i*D + 1];
    if constexpr (D==3)
    {
      file << ' ' << particles.velocities[i*D + 1] <<"\n";
    }
    else
    {
      file <<" 0.0\n";
    }
  }
  file << "        </DataArray>\n";

  file << "        <DataArray type=\"Float64\" Name=\"forces\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int i=0; i<particles.nb_particles; ++i)
  {
    file << "          " << std::setprecision(17) << particles.forces[i*D + 0] << ' ' << particles.forces[i*D + 1];
    if constexpr (D==3)
    {
      file << ' ' << particles.forces[i*D + 1] <<"\n";
    }
    else
    {
      file <<" 0.0\n";
    }
  }
  file << "        </DataArray>\n";
  file << "      </PointData>\n";

  file << "      <Verts>\n";
  file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (std::size_t i = 0; i < particles.nb_particles; ++i) {
    if (i > 0) {
      file << ' ';
    }
    file << i;
  }
  file << "\n        </DataArray>\n";

  file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i = 1; i <= particles.nb_particles; ++i) {
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

#endif