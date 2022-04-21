#pragma once

#include "parameters.hpp"
#include "vector.hpp"
#include "system.hpp"

class Output {
 public:
     std::ofstream energy;
     std::ofstream energy_fluctuation;
     std::ofstream velocities;
     std::ofstream positions;

     std::ofstream temperature;
     std::ofstream mean_free_time;
     std::ofstream parameters;

     std::fstream readme;

     std::string energy_name;
     std::string energy_fluctuation_name;
     std::string velocities_name;
     std::string counter_name;
     std::string temperature_name;
     std::string readme_name;
     std::string mean_free_time_name;
     std::string positions_name;
     std::string parameters_name;

     Output(std::ifstream &file_names, const Parameters &param);

     void finalize();

     void writeResults(const System &sys);

     // Writes position and velocities to file
     void writePositions(const System &sys);

     void writeVelocities(const System &sys);

     void writeEnergyFluctuation(const System &sys);
};
