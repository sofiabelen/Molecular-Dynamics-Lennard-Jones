#pragma once

struct Parameters {
    int n_sim;
    int n_step;
    int n_stable;
    int n_part;
    int dim;

    double dt;
    double density;
    double velocity_mean;
    double r_cut2;
};
