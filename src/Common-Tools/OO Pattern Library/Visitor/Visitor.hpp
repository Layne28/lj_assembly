#include <tuple>

class Particle {
public:
    std::pair<double,double> calculateForce();
    virtual void applySolver(Solver&);
    virtual Particle clone();
private:
    double x;
    double y;
};

class ParticleA : Particle{

};

class ParticleB : Particle{

};




class Solver{
    virtual void update(ParticleA*);
    virtual void update(ParticleB*);


private:
    double dt;
};