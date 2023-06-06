// This is a class used to keep track of clusters composed either of monomers or more clusters
// Think of a percolating cluster. each percolating cluster cluster is comprised of 

// Consider a system of sticky active browinan paricels, if there are enough contantacts a group of particels can be considerd a cluster
// and will calculate their move togheter to save time. 
// clusters might be formed of 

// better yet, think of an adaptive grid, made of several smaller grids
// based on some criteria it makes sense to decrtize space with a high reslotion
// but low resolution is preferable when possible for speed

// a way to implement this would be a set of grids 


// hmm does this want to be a bunch of points/particels with differnt movement rules?
// need something is multiple sclaes



// OOOh OOh maybe this is a good way of defining odd boundary conditions, by gluing togther rectangles and such? 

// look you've got something made up of some other stuff;
// maybe a bunch of vessicels?


// OOOOH A Ciruct is great for this!

#include<tuple>
#include<vector>
class CircuitComponent {
public:
    virtual std::pair<double,double> impedance();
    virtual double voltage();
    virtual double current();

};
class Capacitor: CircuitComponent {
public:
    virtual std::pair<double,double> impedance(){return}
private:
    double _C

}
class CompositeCircuitComponent : CircuitComponent{
public:
    virtual std::pair<double,double> impedance();
    virtual double voltage();
    virtual double current();
private:
    vector<CircuitComponent*> _components;
}






// class Lattice {
//     Lattice();
// private:
//     double _xbounds[2];
//     double _ybounds[2];
// };

// class SquareLattice: public Lattice {
// public:
//     SquareLattice();

// private:
//     double _xbounds[2];
//     double _ybounds[2];
//     double _h;
// };

// class CompositeLattice : public Lattice{



// };






// class Particle {



// private:
//     double _x;
//     double _y;
// };

