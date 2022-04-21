#pragma once

class Vector {
    public:

        int dim;

        std::vector<double> coord;

        // 3-dim constructor, all alements are set to 0
        Vector();

        // n-dim constructor, all alements are set to 0
        Vector(const int &dim);

        // n-dim constructor, all alements are set to val
        Vector(const int &dim, const double &val);

        // 3-dim constructor
        Vector(const int &x, const int &y, const int& z);

        Vector operator +(const Vector &b);
        Vector operator -(const Vector &b);
        Vector operator *(const Vector &b);
        Vector operator *(const double &scalar);
        Vector operator /(const double &scalar);
        Vector& operator *=(const double &scalar);
        Vector& operator /=(const double &scalar);
        Vector& operator +=(const Vector &b);
        Vector& operator -=(const Vector &b);

        // Returns length squared
        double len_sq();

        // Apply periodic boundary conditions
        void wrap(const double &size);

        // Print coordinates to standar output
        void show();
};

// class System {
//  public:
//      /* n_part        : number of particles
//       * n_step        : number of steps
//       * iteration     : current iteration number
//       * size          : size of unit cell in LJ units
//       * halfsize      : half of size
//       * temperature   : temperature in LJ units = 2 * K / dim
//       * kinetic       : kinetic energy K = \sum 1/2 * v**2
//       * potential     : potential energy P = \sum P_i
//       * energy        : total energy E = K + P
//       * dt            : time step
//       * dt2           : time step squared
//       * mft           : mean free time
//       * mean_energy   : mean energy
//       * velocity_mean : mean velocty
//       * standard_dev_energy: energy standard deviation
//       * e_cut         : correction to the potential
//       * r_cut2        : squared cut radius
//       */
// 
//      int n_part;
//      int n_step;
//      int dim;
//      int iteration;
// 
//      double size;
//      double halfsize;
//      double temperature;
//      double kinetic;
//      double potential;
//      double density;
//      double energy;
//      double dt;
//      double dt2;
//      double mft;
//      double mean_energy;
//      double velocity_mean;
//      double standard_dev_energy;
//      double e_cut;
//      double r_cut2;
// 
//      std::vector<Vector> pos;
//      std::vector<Vector> pos_unwrap;
//      std::vector<Vector> vel;
//      std::vector<Vector> acc;
//      std::vector<double> total_energy;
// 
//      void initPositions(
//              int &u,
//              int k,
//              std::vector<int> d,
//              const int &m,
//              const double &ds);
// 
//      void initPositions3d(const int &m, const double &ds);
// 
//      System(
//              const int &n_part,
//              const int &n_step,
//              const int &dim,
//              const double &dt,
//              const double &density,
//              const double &velocity_mean,
//              const double &r_cut2);
// 
//      void calculate_forces();
//      void calculations();
//      void calculateEnergyFluctuation();
//      void show();
// };
