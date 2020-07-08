#include <stdio.h>
#include <ctime>
#include <climits>
#include <cfloat>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>

class Vector {
    public:
        int dim;
        std::vector<double> coord;

        Vector() {
            this->dim = 3;
            coord.resize(this->dim, 0);
        }
        Vector(const int &dim) {
            this->dim = dim;
            coord.resize(dim, 0);
        }

        Vector(const int &dim, const double &val) {
            this->dim = dim;
            coord.resize(dim, val);
        }

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

        void wrap(const double &size) {
            for (int i = 0; i < dim; i++) {
                if (coord[i] < 0) {
                    coord[i] += size;
                } else if (coord[i] >= size) {
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
     int n_part, dim, counter;
     double size, temperature, kinetic, potential, density;
     double energy, dt, mft;
     std::vector<Vector> pos, pos_unwrap, vel, acc, lattice;

     void createLattice(int j, Vector cur) {
         if (j == dim) {
             lattice.push_back(cur);
         } else {
             for (int k = -1; k <= 1; k++) {
                 cur.coord[j] = k * size;
                 createLattice(j + 1, cur);
             }
         }
     }

     void initPositions(int &u, int k, std::vector<int> d,
             const int &m, const double &ds) {
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

     System(const int &n_part, const int &dim, const double &dt,
             const double &density, const double &vel_range) {
         counter = 0;

         this->n_part = n_part;
         this->dim = dim;
         this->density = density;
         this->dt = dt;
         size = pow(n_part / density, 1.0/dim);

         pos.resize(n_part, Vector(dim));
         vel.resize(n_part, Vector(dim));
         acc.resize(n_part, Vector(dim, 0));
         pos_unwrap.resize(n_part, Vector(dim));

         std::default_random_engine generator(time(0));
         std::uniform_real_distribution<double> pos_distribution(0, size);
         std::uniform_real_distribution<double> vel_distribution(-vel_range,
                 vel_range);

         for (int i = 0; i < n_part; i++) {
             for(int j = 0; j < dim; j++) {
                vel[i].coord[j] = vel_distribution(generator);
             }
         }

         int m = static_cast<int>(pow(static_cast<double>(n_part),
                     1.0 / dim));

         // if (static_cast<int>(pow(static_cast<double>(m), dim))
         //     < n_part) {
         //     m++;
         // }

         std::vector<int> d(dim, 0);
         printf("m:%d\n", m);
         int u = 0;
         initPositions(u, 0, d, m, size / static_cast<double>(m + 1));

         pos_unwrap = pos;

         createLattice(0, Vector(dim, 0));
     }

     void calculate_forces() {
// First Step Leapfrog: v(t + 0.5*dt) = v(t) + (0.5*dt)a(t)
//                      r(t + dt) = r(t) + dt * v(t + 0.5*dt)
         for (int i = 0; i < n_part; i++) {
             vel[i] += acc[i] * (0.5 * dt);
             pos[i] += vel[i] * dt;
             pos_unwrap[i] += vel[i] * dt;
             pos[i].wrap(size);
             acc[i].coord.assign(dim, 0);
         }

         potential = 0;
         for (int i = 0; i < n_part; i++) {
            for (int j = i + 1; j < n_part; j++) {
                Vector min_dr;
                double min_dist = DBL_MAX;
                for (int k = 0; k < lattice.size(); k++) {
                    Vector r2 = pos[j] + lattice[k];
                    Vector dr = pos[i] - r2;
                    double dist_ij = dr.len_sq();

                    if (min_dist > dist_ij) {
                        min_dist = dist_ij;
                        min_dr = dr;
                    }
                }

                // (r_ij)^(-2)
                double rri = 1.0 / min_dist;

                // (r_ij)^(-6)
                double rri3 = rri * rri * rri;

                // 48 * (r_ij^(-13) - 0.5*r_ij^(-7))
                double force_ij = 48.0 * rri3 * (rri3 - 0.5) *rri;
                acc[i] += min_dr * force_ij;
                acc[j] -= min_dr * force_ij;

                potential += 4.0 * rri3 * (rri3 - 1.0) + 1.0;
            }
         }

// Second Step Leapfrog: v(t + dt) = v(t + 0.5*dt) + (0.5*dt)a(t + dt)
         for (int i = 0; i < n_part; i ++) {
             vel[i] += acc[i] * (0.5 * dt);
         }

         counter++;
     }

     void calculations() {
         double v_sq_sum = 0;

         for (int i = 0; i < n_part; i++) {
             v_sq_sum += vel[i].len_sq();
         }
         // printf("%f\n", v_sq_sum);

        // Temperature = (n * dim)^(-1) * v_sq_sum;
        temperature = v_sq_sum / static_cast<double>(dim * n_part);
        kinetic = v_sq_sum / (2.0 * static_cast<double>(n_part));
        potential /= static_cast<double>(n_part);
        
        // Mean free time
        mft = 1.0 / (pow(2, 4.0 / 3.0) * pow(kinetic, 1.0 / 3.0) * M_PI * density);
     }

     void show() {
         printf("Number of particles: %d\n", n_part);
         printf("Dimensions: %d\n", dim);
         printf("Size: %.1f\n", size);
         printf("dt: %.1f\n", dt);
         printf("Density: %.1f\n", density);
         printf("Positions:\n");
         for (int i = 0; i < n_part; i++) {
             pos[i].show();
         }
         printf("\nVelocities:\n");
         for (int i = 0; i < n_part; i++) {
             vel[i].show();
         }
         printf("\nPositions unwrapped:\n");
         for (int i = 0; i < n_part; i++) {
             pos_unwrap[i].show();
         }
         printf("\n Lattice (size: %d):\n", lattice.size());
         for (int i = 0; i < lattice.size(); i++) {
             lattice[i].show();
         }
     }
};

struct Parameters {
    int n_sim, n_step, n_stable, n_part, dim;
    double dt, density, vel_range;
};

class Output {
 public:
     std::ofstream energy, velocities, positions, temperature;
     std::ofstream readme, mean_free_time;

     Output(std::ifstream &file_names, const Parameters &param) {
         std::string energy_name, velocities_name , counter_name,
             temperature_name, readme_name, mean_free_time_name,
             positions_name;

         file_names >> energy_name >> velocities_name >> positions_name
             >> counter_name >> temperature_name >> readme_name
             >> mean_free_time_name;

         // std::cout << energy_name << " " << velocities_name << " " <<
         //     counter_name << " " << temperature_name << " " <<
         //     readme_name << " " << mean_free_time_name << " " <<
         //     positions_name << "\n";

         // Keeping track of experiment number
         std::ifstream counter_in(counter_name);
         int exp_count;
         counter_in >> exp_count;
         counter_in.close();

         printf("Experiment number#%d\n", exp_count);

         std::ofstream counter_out(counter_name);
         counter_out << (exp_count + 1);
         counter_out.close();
         // ----------------------------------

         std::string exp_count_str = std::to_string(exp_count);

         energy.open(energy_name + exp_count_str);
         velocities.open(velocities_name + exp_count_str);
         positions.open(positions_name + exp_count_str);
         temperature.open(temperature_name + exp_count_str);
         mean_free_time.open(mean_free_time_name + exp_count_str);

         // double n_part_real = static_cast<double>(param.n_part);
         // velocities << n_part << " " << n_
         // positions << 
         energy << "time kinetic potential\n";
         temperature << "time temp\n";

         std::ofstream readme(readme_name);
         readme << "\nExperiment #" << exp_count << "\n" << param.n_sim
             << "\n" << param.n_step << "\n" << param.n_stable
             << "\n" << param.n_part << "\n" << param.dim << "\n"
             << param.dt << "\n" << param.density << "\n"
             << param.vel_range << "\n\n";
         readme.close();
     }

     void writeResults(const System &sys) {
         energy << sys.counter * sys.dt << " " << sys.kinetic
             << " " << sys.potential << "\n";

         temperature << sys.counter * sys.dt << " "
             << sys.temperature << "\n";

         mean_free_time << sys.counter * sys.dt << " "
             << sys.mft << "\n";
     }

     void writePositions(const System &sys) {
         positions << sys.counter * sys.dt <<"\n";

         for (int i = 0; i < sys.n_part; i++) {
             for (int j = 0; j < sys.dim; j++) {
                positions << sys.pos_unwrap[i].coord[j] << " ";
             }
             positions << "\n";
         }
     }

     void writeVelocities(const System &sys) {
         velocities << sys.counter * sys.dt <<"\n";

         for (int i = 0; i < sys.n_part; i++) {
             for (int j = 0; j < sys.dim; j++) {
                velocities << sys.vel[i].coord[j] << " ";
             }
             velocities << "\n";
         }
     }
};

int main(int argc, char** argv) {
    Parameters param;

    std::ifstream input("Parameters"), file_names("FileNames");
    input >> param.n_sim >> param.n_step >> param.n_stable
        >> param.n_part >> param.dim >> param.dt >> param.density
        >> param.vel_range;

    Output output(file_names, param);

    for (int k = 0; k < param.n_sim; k++) {
        printf("Simulation #%d\n", k);
        System sys(param.n_part, param.dim, param.dt,
                param.density, param.vel_range);
        // sys.show();

        for (int i = 0; i < param.n_step; i++) {
            sys.calculate_forces();
            sys.calculations();
            output.writeResults(sys);

            if (i == param.n_stable) {
                for (int j = 0; j < param.n_part; j++) {
                    sys.pos_unwrap[j] = sys.pos[j];
                }
            }
            if (i >= param.n_stable) { 
                output.writePositions(sys);
            }
        }
        output.writeVelocities(sys);
    }
    return 0;
}
