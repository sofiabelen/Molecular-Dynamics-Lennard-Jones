#include <stdio.h>
#include <ctime>
#include <iostream>
#include <climits>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <ios>

#include "vector.hpp"
#include "system.hpp"
#include "parameters.hpp"
#include "output.hpp"
#include "gui.hpp"

int main(int argc, char** argv) {

    Parameters param;

    std::ifstream input("../Parameters"), file_names("../FileNames");

    input >> param.n_sim
        >> param.n_step
        >> param.n_stable
        >> param.n_part
        >> param.dim
        >> param.dt
        >> param.density
        >> param.velocity_mean
        >> param.r_cut2;

    Output output(file_names, param);

    int window_length = 600;
    int border        = 50;
    int wait_time     = 20;

    for (int k = 0; k < param.n_sim; k++) {

        printf("Simulation #%d\n", k);

        System sys(
                param.n_part,
                param.n_step,
                param.dim,
                param.dt,
                param.density,
                param.velocity_mean,
                param.r_cut2);

        GUI gui(window_length, border, wait_time, sys);

        for (int i = 0; i < param.n_step; i++) {

            if (i == param.n_stable) {
                for (int j = 0; j < param.n_part; j++) {
                    sys.pos_unwrap[j] = sys.pos[j];
                }
            }
            if (i >= param.n_stable) { 
                if (k == 0) {
                    output.writePositions(sys);
                }
                output.writeVelocities(sys);
            }

            gui.redraw(sys);

            sys.calculate_forces();
            sys.calculations();
            output.writeResults(sys);
        }

        sys.calculateEnergyFluctuation();
        output.writeEnergyFluctuation(sys);
        sys.show();

        gui.close();
    }

    output.finalize();
    return 0;
}
