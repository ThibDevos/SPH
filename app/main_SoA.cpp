
#include <particle.h>
#include <initialisation.h>
#include <annexes_functions.h>
#include <kernel_functions.h>
#include <sph_functions.h>
#include <morton_AoS.h>
#include <morton_SoA.h>
#include <write.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>




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
  if(neighbouring!="morton" && neighbouring !="naive")
    std::runtime_error("Choose between morton and naive");
  
  sph::configure(particle_count);
  
  std::cout << "Running simulation with " << particle_count << " target particles"
  << " (actual grid: " << sph::constants().actual_particle_count << ")"
  << std::endl;
  
  file_name_vtp = neighbouring + "_SoA_"+ std::to_string(sph::constants().actual_particle_count);
  file_name_time = file_name_vtp + ".dat";
  std::ofstream f_time;
  f_time.open(file_name_time);

  Particles<2> particles(sph::constants().actual_particle_count);

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
    std::cout<<"filling neighbours"<<std::endl;
    for(int i = 0; i<particles.nb_particles; ++i)
    {
      std::cout<<i<<std::endl;
      int j=0;
      particles.neighbours[i].resize(particles.nb_particles);
      std::iota(particles.neighbours[i].begin(), particles.neighbours[i].end(), 0);
    }
  }

  double T = 2;
  double t = 0.0;
  int count = 0;
  int total_steps = 2000;
  int count_files = 0;
  int save_freq = 200;

  cubic_spline<2> kernel_init(sph::constants().H);

  while (count < total_steps)
  {
    std::cout << count << "/" << total_steps << std::endl;

    if(neighbouring=="morton")
    {
      begin = std::chrono::steady_clock::now();
        Morton_SoA<2> M(particles, sph::constants().H);
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
    
    f_time << count<<" "<<time_initialisation<<" "<<time_computation<<" "<<time_initialisation+time_computation<<"\n";

    update(particles);
    check_boundaries(particles, dom_bot_left, dom_top_right);
    // if (count % save_freq == 0)
    // {

    //    const std::filesystem::path output_path =
    //        std::filesystem::path(__FILE__).parent_path().parent_path() / "bin" / (file_name_vtp + "_"+ std::to_string(count_files) + ".vtp");
    //        count_files++;
    //        write_particles_vtu(particles, output_path.string());
    //  }
     ++count;
    t += sph::constants().DT;
  }
  f_time.flush();
}
