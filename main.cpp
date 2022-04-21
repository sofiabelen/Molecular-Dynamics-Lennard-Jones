#include <stdio.h>
#include <ctime>
#include <iostream>
#include <climits>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <ios>
#include <eggx.h>

class Vector {
    public:

        int dim;

        std::vector<double> coord;

        // 3-dim constructor, all alements are set to 0
        Vector() {

            this->dim = 3;

            coord.resize(this->dim, 0);
        }

        // n-dim constructor, all alements are set to 0
        Vector(const int &dim) {

            this->dim = dim;

            coord.resize(dim, 0);
        }

        // n-dim constructor, all alements are set to val
        Vector(const int &dim, const double &val) {

            this->dim = dim;

            coord.resize(dim, val);
        }

        // 3-dim constructor
        Vector(const int &x, const int &y, const int& z) {

            dim = 3;

            coord.resize(dim);

            coord[0] = x;
            coord[1] = y;
            coord[2] = z;
        }

        Vector operator +(const Vector &b) {

            Vector sum(this->dim);

            for (int i = 0; i < this->dim; i++) {
                sum.coord[i] = this->coord[i] + b.coord[i];
            }
            return sum;
        }

        Vector operator -(const Vector &b) {

            Vector subtr(this->dim);

            for (int i = 0; i < this->dim; i++) {
                subtr.coord[i] = this->coord[i] - b.coord[i];
            }
            return subtr;
        }

        Vector operator *(const Vector &b) {

            Vector dotted(this->dim);

            for (int i = 0; i < this->dim; i++) {
                dotted.coord[i] = this->coord[i] * b.coord[i];
            }
            return dotted;
        }

        Vector operator *(const double &scalar) {

            Vector product(this->dim);

            for (int i = 0; i < this->dim; i++) {
                product.coord[i] = this->coord[i] * scalar;
            }
            return product;
        }

        Vector operator /(const double &scalar) {

            Vector div(this->dim);

            for (int i = 0; i < this->dim; i++) {
                div.coord[i] = this->coord[i] / scalar;
            }
            return div;
        }

        Vector& operator *=(const double &scalar) {

            for (int i = 0; i < this->dim; i++) {
                this->coord[i] *= scalar;
            }
            return *this;
        }

        Vector& operator /=(const double &scalar) {

            for (int i = 0; i < this->dim; i++) {
                this->coord[i] /= scalar;
            }
            return *this;
        }

        Vector& operator +=(const Vector &b) {
            for (int i = 0; i < dim; i++) {
                coord[i] += b.coord[i];
            }
            return *this;
        }

        Vector& operator -=(const Vector &b) {
            for (int i = 0; i < dim; i++) {
                coord[i] -= b.coord[i];
            }
            return *this;
        }

        double len_sq() {
            double sum = 0;
            for (int i = 0; i < dim; i++) {
                sum += coord[i] * coord[i];
            }
            return sum;
        }

        // Apply periodic boundary conditions
        void wrap(const double &size) {

            for (int i = 0; i < dim; i++) {

                if (coord[i] < 0) {
                    coord[i] += size;
                } else if (coord[i] > size) {
                    coord[i] -= size;
                }
            }
        }

        void show() {
            for (int i = 0; i < dim; i++) {
                printf("%.2f ", this->coord[i]);
            }
            printf("\n");
        }
};

class System {
 public:
     /* n_part        : number of particles
      * n_step        : number of steps
      * iteration     : current iteration number
      * size          : size of unit cell in LJ units
      * halfsize      : half of size
      * temperature   : temperature in LJ units = 2 * K / dim
      * kinetic       : kinetic energy K = \sum 1/2 * v**2
      * potential     : potential energy P = \sum P_i
      * energy        : total energy E = K + P
      * dt            : time step
      * dt2           : time step squared
      * mft           : mean free time
      * mean_energy   : mean energy
      * velocity_mean : mean velocty
      * standard_dev_energy: energy standard deviation
      * e_cut         : correction to the potential
      * r_cut2        : squared cut radius
      */

     int n_part;
     int n_step;
     int dim;
     int iteration;

     double size;
     double halfsize;
     double temperature;
     double kinetic;
     double potential;
     double density;
     double energy;
     double dt;
     double dt2;
     double mft;
     double mean_energy;
     double velocity_mean;
     double standard_dev_energy;
     double e_cut;
     double r_cut2;

     std::vector<Vector> pos;
     std::vector<Vector> pos_unwrap;
     std::vector<Vector> vel;
     std::vector<Vector> acc;
     std::vector<double> total_energy;

     void initPositions(
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
    
     void initPositions3d(const int &m, const double &ds) {
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

     System(
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

     void calculate_forces() {

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

     void calculations() {
        temperature = 2.0 * kinetic / static_cast<double>(dim * n_part);
        mean_energy += kinetic + potential;
        total_energy.push_back(kinetic + potential);
        
        // Mean free time
        // mft = 1.0 / (pow(2, 4.0 / 3.0) * pow(kinetic, 1.0 / 3.0) * M_PI * density);
     }

     void calculateEnergyFluctuation() {
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

     void show() {
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
};

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

     Output(std::ifstream &file_names, const Parameters &param) {

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

     void finalize() {
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

     void writeResults(const System &sys) {

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
     void writePositions(const System &sys) {

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

     void writeVelocities(const System &sys) {
         for (int i = 0; i < sys.n_part; i++) {
             for (int j = 0; j < sys.dim; j++) {
                velocities << sys.vel[i].coord[j] << " ";
             }
             velocities << "\n";
         }
     }

     void writeEnergyFluctuation(const System &sys) {
         energy_fluctuation << sys.dt << " " << sys.standard_dev_energy << "\n";
     }
};

class GUI {
    private:
        int win;
        int window_length;
        int border;
        int wait_time;

        double scale;

        std::string bgdark="#212121";
        std::string bglight="#B2DFDB";
        std::string accent="#00BCD4";
        std::string light="#BDBDBD";
        std::string primary="#009688";
        std::string divider="#BDBDBD";
        std::string darkprimary="#212121";

    public:
        GUI(int window_length, int border, int wait_time, System sys) {

            this->window_length = window_length;
            this->border        = border;
            this->wait_time     = wait_time;

            double volume = pow(double(window_length), 3.0);
            double mass   = volume * sys.density;

            scale = pow(mass / double(sys.n_part), 1.0 / 3.0);

            win = gopen(
                    window_length + border * 2,
                    window_length + border * 2);

            layer(win, 0, 1);

            gsetbgcolor(win, bgdark.c_str());
        }

        void redraw(System sys) {
            gclr(win);

            newcolor(win, divider.c_str());

            drawrect(
                    win,
                    double(border),
                    double(border),
                    double(window_length),
                    double(window_length));

            newcolor(win, primary.c_str());

            for(int i=0; i < sys.n_part; i++) {
                // Projection of the particle position onto a 2d plane
                fillcirc(
                        win,
                        sys.pos[i].coord[0] * scale + double(border),
                        sys.pos[i].coord[1] * scale + double(border),
                        5.0,
                        5.0);
            }

            copylayer(win, 1, 0);

            if(wait_time)
                msleep(wait_time);
        }

        void close() {
            // Get keyboard input from user and close the window
            ggetch();
            gclose(win);
        }
};

int main(int argc, char** argv) {

    Parameters param;

    std::ifstream input("Parameters"), file_names("FileNames");

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
