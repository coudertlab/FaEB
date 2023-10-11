#include <iostream>      // printing
#include <string>
#include <chrono>      // timer

#include <local/gridcalc.hpp>
#include <local/clustering.hpp>


int main(int argc, char* argv[]) {
  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
  std::string structure_file = argv[1];
  std::string forcefield_path = argv[2];
  double temperature = std::stod(argv[3]);
  double cutoff = std::stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  std::string element_guest_str = argv[5];
  double approx_spacing = std::stod(argv[6]);
  double energy_threshold = 40;
  double access_coeff = 0.8;
  if (argv[8]==NULL) {energy_threshold = std::stod(argv[7]);}
  if (argv[9]==NULL) {energy_threshold = std::stod(argv[7]);access_coeff = std::stod(argv[8]);}

  // Error catch
  if ( temperature < 0 ) {throw std::invalid_argument( "Received negative value for the Temperature" );}
  if ( energy_threshold < 0 ) {throw std::invalid_argument( "Received negative value for the Energy Threshold" );}
  if ( access_coeff > 1 || access_coeff < 0 ) {throw std::invalid_argument( "Accessibility Coefficient out of range (Read the purpose of this coeff)" );}

  // key constants
  double const R = 8.31446261815324e-3; // kJ/mol/K
  double const N_A = 6.02214076e23;    // part/mol

  // Input
  double molar_mass = 0;
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  gemmi::Grid<double> grid;

  // read_grid_ccp4(grid, grid_file);

  make_energy_grid_ads(structure_file,forcefield_path,temperature,cutoff,element_guest_str,approx_spacing,
  energy_threshold, access_coeff, molar_mass, boltzmann_energy_lj, sum_exp_energy, grid, true);
  const size_t V = grid.data.size();
  
  // //Breadth first search to get the connected components
  size_t N_label = 0;
  uint8_t* channel_labels = new uint8_t[V]();
  conn_comp_energy_bfs(&channel_labels, vector<bool>(V, false), grid, energy_threshold, V, N_label);

  // Array of the type of channel connectivity in X Y Z but not the dimensionality (BFS to do so)
  // Used to filter out pockets no 
  vector<std::string> channel_dimensions=channel_dim_array<uint8_t>(channel_labels, N_label, grid.nu, grid.nv, grid.nw);
  vector<uint8_t> channels;
  for (uint8_t label=0; label!=N_label; label++) { 
    if (channel_dimensions[label]!="\0") {
      // std::cout << label + 1 << " " << channel_dimensions[label] << std::endl;
      channels.push_back(label+1);
    }
  }
  // std::cout << channels.size() << " channels out of " << N_label << " connected clusters" << std::endl;

  // Vector of channel labels grouped by symmetry
  bool error = false;
  vector < vector<uint8_t> > channel_unique_labels = sym_unique_labels(grid, channel_labels, channels, error, std::min(0.0,energy_threshold));
  // print_unique_labels(channel_unique_labels);

  // Loop over the different energy levels
  // double energy_step = R*temperature/3;
  double energy_step = 0.1; // kJ/mol
  vector<std::pair<double,double>> energy_barriers;
  for (auto labels: channel_unique_labels){
    auto label = labels[0];
    double min_energy = energy_threshold; 
    vector<bool> in_channel(V, true); 
    setup_channel_config(in_channel, channel_labels, label, V, grid.data, min_energy);
    // std::cout << min_energy << std::endl;
    
    double energy_threshold_temp = min_energy + energy_step;
    size_t max_steps = floor((energy_threshold-energy_threshold_temp)/energy_step);
    size_t N_current; size_t N_past=0;
    // calc weight here later (TODO)
    uint8_t* bassin_labels_current = new uint8_t[V]();
    bool first=true;
    for (size_t step=0; step<max_steps+1; step++){
      energy_threshold_temp += energy_step;
      N_current = 0;
      bool merged = false;
      conn_comp_energy_bfs(&bassin_labels_current, in_channel, grid, energy_threshold_temp, V, N_current);
      if (first) {
        vector<std::string> channel_dimensions_temp=channel_dim_array<uint8_t>(bassin_labels_current, N_current, grid.nu, grid.nv, grid.nw);
        for (size_t label=0; label!=N_current; label++){
          if (!channel_dimensions_temp[label].empty()) {energy_barriers.push_back(std::make_pair(min_energy, energy_threshold_temp)); first = false; break;}
        }
      }
      N_past = N_current;
    }
    delete [] bassin_labels_current;
  }

  delete [] channel_labels;
  std::string structure_name = trim(structure_file);
  double Framework_density = molar_mass/(N_A*grid.unit_cell.volume*1e-30); // g/m3
  double enthalpy = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  double henry = sum_exp_energy/(grid.data.size()*R*temperature*Framework_density); // mol/kg/Pa
  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  if (energy_barriers.size()!=0) {
    for (auto energy: energy_barriers){
      // Structure name, Enthalpy (kJ/mol), Henry constant (mol/kg/Pa), Minimum Energy in channel (kJ/mol), Barrier Energy (kJ/mol), Time (s)
      std::cout << structure_name << "," << enthalpy << "," << henry << "," << energy.first << "," << energy.second << "," << elapsed_time_ms*0.001 << std::endl;
    }
  }
  else {std::cout << structure_name << "," << enthalpy << "," << henry << ",nan,nan," << elapsed_time_ms*0.001 << std::endl;}
  if (error) {std::cout << structure_name << ": error" << std::endl;}
}
