//A System consists of a set of Particles along with boundary conditions and control parameters (temperature, external field, etc.)

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include "Particle.hpp"
#include "ParamDict.hpp"
#include "Observer.hpp"

class Observer;

class System
{
private:
    Observer *obs; //use to write output

public:
    //For now we assume NVT ensemble
    int N; //No. of particles (nodes)
    double rho; //density
    double phi; //packing fraction
    double kT; //temperature
    int dim; //# of spatial dimensions
    double dt; //timestep

    double rcut; //cutoff distance for pair potential
    double sigma;
    double epsilon;

    int time; //No. of timesteps taken (can be reset)

    std::vector<Particle> particles; //all the particles in the system
    std::vector<std::vector<int>> image; //which periodic image each particle is in (starts at (0,0,0))
    std::vector<double> edges; //box dimensions
    std::vector<int> is_periodic; //periodic or not in each dimension

    /*** Methods ***/

    //constructor
    //TODO: add default ParamDict
    System(ParamDict &theParams);

    //destructor
    ~System();

    //Make changes to System state
    void apply_pbc(); //apply periodic boundary conditions
    void zero_com(); //subtract center of mass from each particle position
    void set_obs(Observer &anObs);

    //Get (some of these could be made static)
    arma::vec get_com();
    double get_energy();
    double get_lj_potential(double r, double sig, double eps, double rc);
    arma::vec get_force(Particle &p);
    arma::vec get_lj_force(double r, arma::vec rvec, double sig, double eps);
    arma::vec get_disp_vec(Particle &p1, Particle &p2);
    double get_dist(Particle &p1, Particle &p2);
    Observer get_obs();
};

#endif