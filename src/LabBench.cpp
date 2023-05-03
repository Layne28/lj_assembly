#include "LabBench.hpp"

LabBench::LabBench(ParamDict& theParams, gsl_rng*& theGen) : sys(theParams), obs(theParams), solver(sys, theParams, theGen)
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
    std::cout << "Equilibrating..." << std::endl;
    this->run(this->equil_steps, "/equil", this->obs.particles_freq, this->obs.thermo_freq);

    std::cout << "Doing production run..." << std::endl;
    this->run(this->production_steps, "/prod", this->obs.particles_freq, this->obs.thermo_freq);
}
