#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>

#include "output.hpp"

Output::Output(std::ifstream &file_names, const Parameters &param) {

    file_names
        >> energy_name
        >> energy_fluctuation_name
        >> velocities_name
        >> positions_name
        >> counter_name
        >> temperature_name
        >> readme_name
        >> mean_free_time_name
        >> parameters_name;

    // Keeping track of experiment number ------------------------------
    std::ifstream counter_in(counter_name);
    int exp_count;
    counter_in >> exp_count;
    counter_in.close();

    printf("Experiment number #%d\n", exp_count);

    std::ofstream counter_out(counter_name);
    counter_out << (exp_count + 1);
    counter_out.close();
    // -----------------------------------------------------------------

    std::string exp_count_str = std::to_string(exp_count);

    energy.open(energy_name + exp_count_str);
    energy_fluctuation.open(energy_fluctuation_name + exp_count_str);
    velocities.open(velocities_name + exp_count_str);
    positions.open(positions_name + exp_count_str + ".xyz");
    temperature.open(temperature_name + exp_count_str);
    mean_free_time.open(mean_free_time_name + exp_count_str);
    parameters.open(parameters_name + exp_count_str);

    energy << "time kinetic potential\n";
    energy_fluctuation << "dt fluctuation\n";
    temperature << "time temp\n";

    for (int i = 0; i < param.dim; i++) {
        velocities << i << " ";
    }
    velocities << "\n";

   parameters
       << param.n_sim  << "\n"
       << param.n_step        << "\n"
       << param.n_stable      << "\n"
       << param.n_part        << "\n"
       << param.dim           << "\n"
       << param.dt            << "\n"
       << param.density       << "\n"
       << param.velocity_mean << "\n"
       << param.r_cut2;

    readme.open(readme_name, std::ios_base::app);
    readme
        << "\nExperiment #"
        << exp_count           << "\n"
        << param.n_sim         << "\n"
        << param.n_step        << "\n"
        << param.n_stable      << "\n"
        << param.n_part        << "\n"
        << param.dim           << "\n"
        << param.dt            << "\n"
        << param.density       << "\n"
        << param.velocity_mean << "\n"
        << param.r_cut2        << "\n\n";
    readme.close();
}

void Output::finalize() {
    readme.open(readme_name, std::ios_base::app);
    readme << "Experiment finalized.\n";

    readme.close();
    energy.close();
    energy_fluctuation.close();
    velocities.close();
    positions.close();
    temperature.close();
    mean_free_time.close();
}

void Output::writeResults(const System &sys) {

    energy
        << sys.iteration * sys.dt
        << " "
        << sys.kinetic
        << " "
        << sys.potential
        << "\n";

    temperature
        << sys.iteration * sys.dt
        << " "
        << sys.temperature
        << "\n";

    mean_free_time
        << sys.iteration * sys.dt
        << " "
        << sys.mft
        << "\n";
}

// Writes position and velocities to file
void Output::writePositions(const System &sys) {

    // Changed pos_unwrap -> pos
    positions
        << sys.n_part
        << "\n"
        << sys.iteration
        << "\n";

    for (int i = 0; i < sys.n_part; i++) {

        for (int j = 0; j < sys.dim; j++) {
           positions << sys.pos[i].coord[j] << " ";
        }

        for (int j = 0; j < sys.dim; j++) {
           positions << sys.vel[i].coord[j] << " ";
        }

        positions << "\n";
    }
}

void Output::writeVelocities(const System &sys) {
    for (int i = 0; i < sys.n_part; i++) {
        for (int j = 0; j < sys.dim; j++) {
           velocities << sys.vel[i].coord[j] << " ";
        }
        velocities << "\n";
    }
}

void Output::writeEnergyFluctuation(const System &sys) {
    energy_fluctuation << sys.dt << " " << sys.standard_dev_energy << "\n";
}
