#include <stdio.h>
#include <iostream>
#include <vector>
#include <random>

#include "vector.hpp"
#include "system.hpp"

void System::initPositions(
        int &u,
        int k,
        std::vector<int> d,
        const int &m,
        const double &ds) {

    if (u < n_part) {

       for (int i = 0; i < dim; i++) {
            pos[u].coord[i] = ds * static_cast<double>(d[i]);
       }

       u++;

       for (int j = k; j < dim; j++) {
           if(d[j] + 1 < m) {
               d[j]++;
               initPositions(u, j, d, m, ds);
               d[j]--;
           }
       }
    }
}

void System::initPositions3d(const int &m, const double &ds) {
    int u = 0;
    for(int i = 0; i < m; i++) {
       for(int j = 0; j < m; j++) {
           for(int k = 0; k < m; k++) {

               pos[u].coord[0] = ds * (i + 0.5);
               pos[u].coord[1] = ds * (j + 0.5);
               pos[u].coord[2] = ds * (k + 0.5);

               u++;
           }
       }
    }
}

System::System(
        const int &n_part,
        const int &n_step,
        const int &dim,
        const double &dt,
        const double &density,
        const double &velocity_mean,
        const double &r_cut2) {

    iteration = 0;

    this->n_part        = n_part;
    this->n_step        = n_step;
    this->dim           = dim;
    this->density       = density;
    this->dt            = dt;
    this->r_cut2        = r_cut2;
    this->velocity_mean = velocity_mean;
    this->dt2           = dt * dt;

    pos.resize(n_part, Vector(dim));

    size = pow(n_part / density, 1.0/dim);
    halfsize = size / 2.0;
    double rr3 = 1.0 / (r_cut2 * r_cut2 * r_cut2);
    e_cut = 4 * (rr3 * rr3 * rr3 * rr3 - rr3 * rr3);
    std::cout << "ecut: " << e_cut << "\n";

    pos.resize(n_part, Vector(dim));
    vel.resize(n_part, Vector(dim));
    acc.resize(n_part, Vector(dim, 0));
    pos_unwrap.resize(n_part, Vector(dim));

    std::random_device rd;
    std::mt19937 generator(rd());
    std::exponential_distribution<> vel_distribution(velocity_mean);

    for (int i = 0; i < n_part; i++) {
        for(int j = 0; j < dim; j++) {
           vel[i].coord[j] = vel_distribution(generator);
        }
    }

    // Take away any center-of-mass drift
    Vector cmv;

    for (int i = 0; i < n_part; i++) {
        for(int j = 0; j < dim; j++) {
            cmv.coord[j] += vel[i].coord[j];
        }
    }

    for (int i = 0; i < n_part; i++) {
        for(int j = 0; j < dim; j++) {
            vel[i].coord[j] -= cmv.coord[j] / n_part;
        }
    }

    // Initialize particles in a grid ----------------------------------
    int m = static_cast<int>(pow(
                static_cast<double>(n_part),
                1.0 / dim));

    if (static_cast<int>(pow(static_cast<double>(m), dim)) < n_part) {
        m++;
    }

    initPositions3d(m, size / static_cast<double>(m));

    pos_unwrap = pos;
    //------------------------------------------------------------------
}

void System::calculate_forces() {

    // First integration half-step
    for (int i = 0; i < n_part; i++) {

        Vector dr = vel[i] * dt + acc[i] * (0.5 * dt2);

        pos[i]        += dr;
        pos_unwrap[i] += dr;

        vel[i]        += acc[i] * (0.5 * dt);

        pos[i].wrap(size);

        // acceleration = 0
        acc[i].coord.assign(dim, 0);
    }

    potential = 0;
    kinetic   = 0;

    // Forces calculation
    for (int i = 0; i < n_part; i++) {
       for (int j = i + 1; j < n_part; j++) {

           Vector dr = pos[i] - pos[j];

           // Periodic boundary conditions:
           // Apply the minimum image convention
           for (int k = 0; k < dim; k++) {
               if (dr.coord[k] < -halfsize)
                   dr.coord[k] += size;
               else if (dr.coord[k] > halfsize)
                   dr.coord[k] -= size;
           }

           double r2 = dr.len_sq();
           
           if (r2 < r_cut2) {

               double r6i = 1.0 / (r2 * r2 * r2);

               // 48 * (r_ij^(-13) - 0.5*r_ij^(-7))
               double force_ij = 48 * (r6i * r6i - 0.5 * r6i);
               acc[i] += dr * (force_ij / r2);
               acc[j] -= dr * (force_ij / r2);

               potential += 4.0 * (r6i * r6i - r6i) - e_cut;
           }
       }
    }

    // Second integration half-step
    for (int i = 0; i < n_part; i ++) {

        vel[i]  += acc[i] * (0.5 * dt);

        kinetic += vel[i].len_sq();

    }

    kinetic *= 0.5;
    iteration++;
}

void System::calculations() {
   temperature = 2.0 * kinetic / static_cast<double>(dim * n_part);
   mean_energy += kinetic + potential;
   total_energy.push_back(kinetic + potential);
   
   // Mean free time
   // mft = 1.0 / (pow(2, 4.0 / 3.0) * pow(kinetic, 1.0 / 3.0) * M_PI * density);
}

void System::calculateEnergyFluctuation() {
    mean_energy /= n_step;
    standard_dev_energy = 0;

    for(int i = 0; i < n_step; i++) {
        double energy_diff = total_energy[i] - mean_energy;
        standard_dev_energy += energy_diff * energy_diff;
    }

    standard_dev_energy /= n_step - 1;
    standard_dev_energy = pow(standard_dev_energy, 0.5);
    // standard_dev_energy /= mean_energy;
}

void System::show() {
    printf("Number of particles: %d\n", n_part);
    printf("Dimensions: %d\n", dim);
    printf("Size: %.1f\n", size);
    printf("dt: %f\n", dt);
    printf("Density: %.1f\n", density);
    std::cout << "velocity_mean: " << velocity_mean << "\n";
    std::cout << "r_cut2: " << r_cut2 << "\n";
    printf("Mean Total Energy: %.1f\n", mean_energy);
    std::cout << "Energy Standard Deviation: "
        << standard_dev_energy << "\n";
   printf("\n\n");
}
