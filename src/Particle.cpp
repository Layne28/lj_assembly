#include <iostream>
#include "Particle.hpp"

//TODO: write unit tests for object id/counter
int Particle::counter = 0;

Particle::Particle(int dim) {
    pos.zeros(dim);
    vel.zeros(dim);
    id = Particle::counter;
    Particle::counter++;
}

Particle::~Particle() {}


int Particle::get_id() {
    return id;
}

arma::vec Particle::get_pos() {
    return pos;
}

bool Particle::is_equal(Particle &p){
    if (this->id==p.id) {
        return true;
    }
    else {
        return false;
    }
}

std::ostream& operator<<(std::ostream& os, Particle& p) {
    arma::vec pos = p.get_pos();
    return os << "Particle with id: " << p.get_id() << " and Pos: (" << pos(0) << ", " << pos(1) << ", " << pos(2) << ")";
}