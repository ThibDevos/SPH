
#include <particle.h>
#include <initialisation.h>
#include <annexes_functions.h>
#include <kernel_functions.h>
#include <sph_functions.h>
#include <morton.h>

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

  file << "        <DataArray type=\"Float64\" Name=\"leaf\" format=\"ascii\">\n";
  for (const auto& particle : particles) {
    file << "          " << std::setprecision(17) << particle.leaf << '\n';
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

int main(int argc, char** argv)
{
  std::size_t particle_count = 1000;
  std::string neighbouring = "morton";
  std::string file_name_time;
  std::string file_name_vtp;
  if (argc > 1)
  {
    particle_count = static_cast<std::size_t>(std::stoul(argv[1]));
    if (argc >2)
      neighbouring = argv[2];
  }
  if(neighbouring!="morton" || neighbouring !="naive")
    std::runtime_error("Choose between morton and naive");
  
  sph::configure(particle_count);
  
  std::cout << "Running simulation with " << particle_count << " target particles"
  << " (actual grid: " << sph::constants().actual_particle_count << ")"
  << std::endl;
  
  file_name_vtp = neighbouring + "_"+ std::to_string(sph::constants().actual_particle_count);
  file_name_time = file_name_vtp + ".dat";
  std::ofstream f_time;
  f_time.open(file_name_time);

  std::vector<Particle<2>> particles(sph::constants().actual_particle_count);

  std::array<double, 2> dam_bot_left{0.0, 0.0};
  std::array<double, 2> dam_top_right{0.146, 0.292};
  const double dx = sph::constants().dx;

  std::array<double, 2> dom_bot_left{0.0, 0.0};
  std::array<double, 2> dom_top_right{0.584, 0.35};

  init_grid_block(dam_bot_left, dam_top_right, dx, particles);
 
  std::chrono::steady_clock::time_point begin;
  std::chrono::steady_clock::time_point end;
  auto time_initialisation =  std::chrono::duration_cast<std::chrono::microseconds>(begin - begin).count();
  auto time_computation = time_initialisation;



  //naive
  if(neighbouring=="naive")
  {
    for(auto & p : particles)
    {
      p.neighbours.clear();
      for(auto & q : particles)
      {
        p.neighbours.push_back(&q);
      }
    }
  }

  double T = 2;
  double t = 0.0;
  int count = 0;
  int count_files = 0;
  int save_freq = 500;

  cubic_spline<2> kernel_init(sph::constants().H);

  while (t < T)
  {
    std::cout << t << "/" << T << std::endl;

    if(neighbouring=="morton")
    {
      begin = std::chrono::steady_clock::now();
        Morton<2> M(particles, sph::constants().H);
      end = std::chrono::steady_clock::now();
      time_initialisation = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    }
    else
    {
      time_initialisation = 0.;
    }
    
    
    
    begin = std::chrono::steady_clock::now();
      compute_density_pressure(particles);
      compute_forces(particles);
    end = std::chrono::steady_clock::now();
    time_computation = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    
    f_time << t<<" "<<time_initialisation<<" "<<time_computation<<" "<<time_initialisation+time_computation<<"\n";

    update(particles);
    check_boundaries(particles, dom_bot_left, dom_top_right);
    //  if (count % save_freq == 0)
    //  {

    //    const std::filesystem::path output_path =
    //        std::filesystem::path(__FILE__).parent_path().parent_path() / "bin" / (file_name_vtp + "_"+ std::to_string(count_files) + ".vtp");
    //        count_files++;
    //        write_particles_vtu(particles, output_path.string());
    //  }
    //  ++count;
    t += sph::constants().DT;
  }
  f_time.flush();
}
