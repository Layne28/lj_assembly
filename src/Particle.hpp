// A Particle object has a position (x, y, z) and a velocity (vx, vy, vz) in 3D space.

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <array>
#include <vector>
#include <armadillo>

class Particle
{
private:
    static int counter; //Warning -- not thread-safe!!
    int id;
    
public:
    arma::vec pos;
    arma::vec vel;

    //constructor
    Particle(int dim=3);

    //destructor
    ~Particle();

    int get_id();
    arma::vec get_pos();
    bool is_equal(Particle &p);

    //'<<' overload for std::cout
    friend std::ostream& operator<<(std::ostream& os, Particle& p);
};

#endif