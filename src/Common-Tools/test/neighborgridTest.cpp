#include "neighborgrid.h"
#include <boost/test/unit_test.hpp>
#include <random>
BOOST_AUTO_TEST_SUITE(neighborgrid_suite);

struct F {

    F(){
        std::array<double,2> minBnd {0,0};
        std::array<double,2> maxBnd {1,1};
        std::array<bool,2> periodic {true,true};
        double r = 0.1;
        g = NeighborGrid<Atom<2>,2>(minBnd,maxBnd,periodic,r);

    }

    NeighborGrid<Atom<2>,2> g;
};

BOOST_FIXTURE_TEST_CASE(two_particles,F){
    double r[3] {0.5,0.5};
    double r2[3] {0.45,0.45};
    Atom<2> particle1(r);
    Atom<2> particle2(r2);


    g.update_atom(&particle1);
    g.update_atom(&particle2);

    std::vector<Atom<2>*> n = g.get_neighbors(&particle1);

    BOOST_CHECK_EQUAL(1,n.size());
    BOOST_CHECK_EQUAL(&particle2,n[0]);

    

}
BOOST_FIXTURE_TEST_CASE(two_particles_periodic,F){
    double r[3] {0.5,0.98};
    double r2[3] {0.5,0.01};
    Atom<2> particle1(r);
    Atom<2> particle2(r2);


    g.update_atom(&particle1);
    g.update_atom(&particle2);

    std::vector<Atom<2>*> n = g.get_neighbors(&particle1);

    BOOST_CHECK_EQUAL(1,n.size());
    BOOST_CHECK_EQUAL(&particle2,n[0]);

}

BOOST_AUTO_TEST_CASE(class_on_top_of_adam){


    class DiffusiveAtom : public Atom<2>{
    protected:
        double _D = 0.0001;
    public:
        DiffusiveAtom(double *r,double D) : Atom<2>(r),_D(D){};
        void diffuse(double dt){
                std::default_random_engine generator;
                std::normal_distribution<double> distribution(0,sqrt(_D)*dt);

            double dr[2];
            dr[0] = distribution(generator);
            dr[1] = distribution(generator);
            move(dr);
        }
    };


    std::array<double,2> minBnd {0,0};
    std::array<double,2> maxBnd {1,1};
    std::array<bool,2> periodic {true,true};
    double r = 0.1;
    NeighborGrid<DiffusiveAtom,2> g1(minBnd,maxBnd,periodic,r);


    double r1[3] {0.5,0.5};
    double r2[3] {0.5,0.5};
    DiffusiveAtom particle1(r1,1e-3);
    DiffusiveAtom particle2(r2,1e-3);
    particle1.diffuse(1e-6);
    particle2.diffuse(1e-6);


    g1.update_atom(&particle1);
    g1.update_atom(&particle2);

    std::vector<DiffusiveAtom *> n = g1.get_neighbors(&particle1);

    BOOST_CHECK_EQUAL(1,n.size());
    if (n.size()==1) BOOST_CHECK_EQUAL(&particle2,n[0]);
    ((DiffusiveAtom *)n[0])->diffuse(1e-3);

    

}


BOOST_AUTO_TEST_SUITE_END();

