#ifndef neighborgrid_h
#define neighborgrid_h
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <vector>
#include <stdexcept>
#include <cmath>
template <int D>
class Atom {
protected:
    std::array<double,D> _position;
    std::array<double,D> _oldPosition = {0};

public:
    Atom(){for (int i=0; i<D; i++) {_position[i]=0;_oldPosition[i]=0;}}
    Atom(double* r) {for (int i=0; i<D; i++) _position[i]=r[i];}
    

    double* get_position(){return _position.data();}
    double* get_old_position(){return _oldPosition.data();}

    void move(double *dr) {
        for (int i = 0; i<D; i++){
            _oldPosition[i]=_position[i];
            _position[i]+=dr[i];
        }
    }

};

template <class C,int D>
class NeighborGrid {
protected: 
    // provided by user
    std::array<double,D> _minBoundary;
    std::array<double,D> _maxBoundary;
    std::array<bool,D> _periodic;
    double _interactionLength;
    
    
    std::unordered_map<long, std::unordered_set<C*>> _grid;    
    double _boxSize;
    std::array<int,D> _numBoxes;
    


    // 2 options here , pass back a vector of neighbors or keep a scratch vector on hand to use.
    // lets start with the vector and if we need to speed things up we can use the a scratch space

    void initialize();
    std::array<int,D> position_to_box_number(double *);
    long box_number_to_index(std::array<int,D>);
    long position_to_index(double *);
    double distance(C *particle1, C *particle2);
    void recurse_into_neighbor_grid(std::array<int,D> position, std::array<int,D> neighbor_position,int dimension,std::vector<C*> &neighbors,C* particle);

public:
    NeighborGrid(){};
    NeighborGrid(std::array<double,D> min,std::array<double,D> max,std::array<bool,D> p,double r):_minBoundary(min),_maxBoundary(max),_periodic(p), _interactionLength(r) {initialize();}
    
    std::vector<C*> get_neighbors(C *);
    
    void update_atom(C*);
    void print();
};

template <class C,int D>
inline void NeighborGrid<C,D>::initialize()
{
    //set up box sizes
    _boxSize = _interactionLength/2;
    for (int i = 0; i<D; i++) {
        _numBoxes[i] = floor((_maxBoundary[i]-_minBoundary[i])/_boxSize);
    }
    
}

template <class C,int D>
std::array<int,D> NeighborGrid<C,D>::position_to_box_number(double * position)
{
   std::array<int,D> index;
    for (int i = 0; i<D; i++){
        if (position[i] < _minBoundary[i] || position[i] >=_maxBoundary[i]){
            throw std::domain_error("Atom out of bounds");
        }
        index[i] = (position[i]-_minBoundary[i])/(_maxBoundary[i]-_minBoundary[i])*_numBoxes[i];
    }
    return index;
}
template <class C,int D>
inline long NeighborGrid<C,D>::position_to_index(double *position)
{
    return box_number_to_index(position_to_box_number(position));
}

template <class C,int D>
inline long NeighborGrid<C,D>::box_number_to_index(std::array<int, D> position)
{
   long index = 0;
   long mult = 1;

   
   for (int i = 0; i < D; i++)
   {
        index += position[i]*mult;
        mult*=_numBoxes[i];
   }
    return index;
}
template <class C,int D>
inline double NeighborGrid<C,D>::distance(C* particle1, C* particle2)
{
    double sum = 0;
    double * pos1 = particle1->get_position();
    double * pos2 = particle2->get_position();
    for (int i = 0; i<D; i++){
        double dist = (pos2[i]-pos1[i]);
        if (_periodic[i]){
            double left = dist - (_maxBoundary[i]-_minBoundary[i]);
            double right = dist + (_maxBoundary[i]-_minBoundary[i]);
            sum+=fmin(left*left,fmin(right*right,dist*dist));
        } else {
            sum+=dist*dist;
        }
    }
    return sqrt(sum);
}


template <class C,int D>
inline void NeighborGrid<C,D>::recurse_into_neighbor_grid(std::array<int,D> position, std::array<int,D> neighbor_position, int dimension, std::vector<C *> &neighbors,C* particle)
{
    for (int i = -2; i<3; i++){
        neighbor_position[dimension]=position[dimension]+i;

        // if we're on the min boundary
        if (neighbor_position[dimension]<0) {
            // if we're not periodic in this dimension continue
            if (!_periodic[dimension]) continue;
            // add the number of boxes to wrap around
            
            neighbor_position[dimension] += _numBoxes[dimension];
        }

        // if we're on the max boundary

        if (neighbor_position[dimension]>=_numBoxes[dimension]) {
            // if we're not periodic in this dimension continue
            if (!_periodic[dimension]) continue;
            // subtract the number of boxes to wrap around
            neighbor_position[dimension]-=_numBoxes[dimension];
        }

        // if we in the lowest dimension then get the particles in this 
        // box and compare their distance, adding them to our neighbor list if 
        // they are close enough
        if (dimension==0) {
            if (_grid.find(box_number_to_index(neighbor_position))!=_grid.end()){
                // if this grid is populated
                std::unordered_set<C *> boxSet = _grid[box_number_to_index(neighbor_position)];
                // loop over particles in it
                for (auto it = boxSet.begin(); it!=boxSet.end(); it++){
                    // if they're in the interaction distance
                    if (*it!=particle && distance(*it,particle)<_interactionLength){
                        // add them
                        neighbors.push_back(*it);
                    }
                }
            }
            
        } else {
            // if we're not in the lowest dimension keep going
            recurse_into_neighbor_grid(position,neighbor_position,dimension-1,neighbors,particle);

        }

    }
}

template <class C,int D>
inline std::vector<C *> NeighborGrid<C,D>::get_neighbors(C * particle)
{
    std::vector<C *> neighbors;
    // first get the position
    std::array<int,D> position = position_to_box_number(particle->get_position());
    std::array<int,D> neighbor_position = position;
    recurse_into_neighbor_grid(position,neighbor_position,D-1,neighbors,particle);
    
    return neighbors;
}

template <class C,int D>
inline void NeighborGrid<C,D>::update_atom(C * particle)
{
    // try to remove this atom from where it was
    // if the box exists
    if (_grid.find(position_to_index(particle->get_old_position()))!=_grid.end()){

        // get the box
        std::unordered_set<C*> &oldBoxSet = _grid[position_to_index(particle->get_old_position())];
        // find the particle in it
        auto particleReference = oldBoxSet.find(particle);
        if (particleReference!=oldBoxSet.end()){ // if we found it
            oldBoxSet.erase(*particleReference);
            // remove this box if its empty
            if (oldBoxSet.size()==0){ 
                _grid.erase(position_to_index(particle->get_old_position()));//questionable to do this here
            }
        }
    }

    // add it to where it should be
    // if we don't have a box there yet, make one
    if (_grid.find(position_to_index(particle->get_position()))==_grid.end()){
        // printf("Adding it to a new box %ld\n",position_to_index(particle->get_position())); 
        
        _grid[position_to_index(particle->get_position())] = std::unordered_set<C*>{particle};

    }
    else
    {
        std::unordered_set<C*> &boxSet = _grid[position_to_index(particle->get_position())];
        boxSet.insert(particle);
    }
}
template <class C,int D>
inline void NeighborGrid<C,D>::print() {
    for (auto it = _grid.begin(); it!=_grid.end(); it++){
        for (auto it2 = (*it).second.begin(); it2 != (*it).second.end(); it2++){
            std::array<double,D> pos = (*it2)->get_position();
            std::array<int,D> index = position_to_box_number(pos);
            if (D==2) {
                printf("Particle %p at [%f,%f] which is index [%d,%d]\n",*it2,pos[0],pos[1],index[0],index[1]);
            } else if (D==3) {
                printf("Particle %p at [%f,%f,%f] which is index [%d,%d,%d]\n",*it2,pos[0],pos[1],pos[2],index[0],index[1],index[2]);
            }

        }
    }
}


#endif