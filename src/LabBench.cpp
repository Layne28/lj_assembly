#include "LabBench.hpp"

LabBench::LabBench(ParamDict& theParams, gsl_rng*& theGen) : sys(theParams, theGen), obs(theParams), solver(sys, theParams, theGen)
{

    params = theParams;
    sys.set_obs(obs);

    if(theParams.is_key("equil_steps")) equil_steps = std::stoi(theParams.get_value("equil_steps"));
    if(theParams.is_key("production_steps")) production_steps = std::stoi(theParams.get_value("production_steps"));
    if(theParams.is_key("info_freq")) info_freq = std::stoi(theParams.get_value("info_freq"));
    if(theParams.is_key("experiment")) experiment = theParams.get_value("experiment");
}

LabBench::~LabBench() {}

void LabBench::run(int nstps, std::string subdir, int net_freq, int therm_freq)
{
    if (nstps==-1) nstps = this->production_steps;
    if (net_freq==-1) net_freq = this->obs.particles_freq;
    if (therm_freq==-1) therm_freq = this->obs.thermo_freq;

    if (obs.do_h5md) {
        obs.open_h5md(sys, subdir);
    }

    for (int i=0; i<nstps; i++) {
        //Record data
        if (i%info_freq==0) std::cout << "step " << i << std::endl;

        if (obs.do_h5md) {
            if (i%net_freq==0) {
                obs.dump_h5md(sys, subdir);
            }
        }
        //Advance dynamics
        solver.update(sys);
    }
}

void LabBench::do_experiment(std::string expt)
{
    if (expt=="standard") {
        std::cout << "Running standard experiment." << std::endl;
        this->run_standard_experiment();
    } 
    else {
        std::cout << "This experiment has not been designed yet.\n" << std::endl;
        exit(0);
    }
}

void LabBench::run_standard_experiment()
{
    std::cout << "Minimizing energy..." << std::endl;
    sys.minimize_energy();

    /*
    std::cout << "Testing grid" << std::endl;
    
    std::array<double,2> min = {-0.5*sys.edges[0], -0.5*sys.edges[1]};
    std::array<double,2> max = {0.5*sys.edges[0], 0.5*sys.edges[1]};
    std::array<bool,2> per = {true, true};
    NeighborGrid<Particle, 2> *grid = new NeighborGrid<Particle,2>(min, max, per, 3.5);
    //grid->update_atom(&sys.particles[0]);
    //grid->update_atom(&sys.particles[1]);
    //grid->update_atom(&sys.particles[2]);
    //auto p = grid->get_neighbor_iterator(&sys.particles[0]);
    //std::cout << *p << std::endl;
    Particle p1;
    p1.pos[0]=0.0;
    p1.pos[1]=0.0;
    Particle p2;
    p2.pos[0] = 1.0;
    p2.pos[1] = 0.0;
    grid->update_atom(&p1);
    grid->update_atom(&p2);
    for (auto p = grid->get_neighbor_iterator(&p1); !p.isDone(); p++){ 
        std::cout << *p << std::endl;
    }
    std::cout << "passed" << std::endl;
    */

    if (sys.particle_protocol=="lattice") {
        std::cout << "Equilibrating..." << std::endl;
        double kT0 = sys.kT;
        double va0 = solver.va;
        double dt0 = sys.dt;
        std::string pot = sys.potential_type;
        
        //Set "high" temperature and no active noise
        //to get random initial configuration
        sys.kT = 0.5;
        solver.va = 0.0;
        sys.dt = 0.00025;
        sys.potential_type = "wca";
        std::cout << "kT: " << sys.kT << std::endl;;
        std::cout << "va: " << solver.va << std::endl;
        std::cout << "potential: " << sys.potential_type << std::endl;
        std::cout << "dt: " << sys.dt << std::endl;;
        this->run(this->equil_steps, "/equil", this->obs.particles_freq, this->obs.thermo_freq);

        //reset parameters for production
        sys.kT = kT0;
        solver.va = va0;
        sys.dt = dt0;
        sys.potential_type = pot;
    }
    std::cout << "Doing production run..." << std::endl;
    std::cout << "kT: " << sys.kT << std::endl;;
    std::cout << "va: " << solver.va << std::endl;
    std::cout << "potential: " << sys.potential_type << std::endl;
    this->run(this->production_steps, "/prod", this->obs.particles_freq, this->obs.thermo_freq);
}
