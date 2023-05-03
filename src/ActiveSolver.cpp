#include "ActiveSolver.hpp"

ActiveSolver::ActiveSolver(System &theSys, ParamDict &theParams, gsl_rng *&the_rg) : Solver(theParams, the_rg)
{
    if(theParams.is_key("va")) va = std::stod(theParams.get_value("va"));
    //Modify parameters to make particles and active noise generator consistent 
    theParams.add_entry("D", theParams.get_value("tau")); //TODO: figure out what the hell this is doing
    double dee_x = theSys.edges[0]/int(std::stoi(theParams.get_value("nx")));
    //Set precision of dx
    std::ostringstream out;
    out.precision(15);
    out << std::fixed << dee_x;
    theParams.add_entry("dx", out.str());
    anGen = new Generator(theParams, the_rg);

    std::cout << anGen->ny << std::endl;

    //Check that active noise gen. params are consistent with particle system
    if(fabs(anGen->Lx-theSys.edges[0])>1e-10) {
        std::cout << "ERROR! System and Active Noise Generator do not have the same size (" << std::fixed << theSys.edges[0] << " vs " << anGen->Lx << ")." << std::endl;
        exit(1);
    }
}

ActiveSolver::~ActiveSolver()
{
    delete anGen;
}

void ActiveSolver::update(System &theSys)
{
    //Get conservative forces
    std::vector<arma::vec> potential_forces(theSys.N);
    for (int i=0; i<theSys.N; i++) {
        arma::vec force = theSys.get_force(theSys.particles[i]);
        potential_forces[i] = force;
    }

    //Get active noise force on each particle
    std::vector<arma::vec> active_forces;
    active_forces = get_active_noise_forces(theSys, *anGen);
    
    for (int i=0; i<theSys.N; i++) {
        arma::vec pos = theSys.particles[i].get_pos();
        for (int k=0; k<theSys.dim; k++) {
            //Euler step
            double incr = potential_forces[i](k)/gamma*dt 
                        + active_forces[i](k)*dt;
            if (theSys.kT>0) incr += sqrt(2*theSys.kT/gamma)*gsl_ran_gaussian(rg, sqrt(dt)); 
            pos(k) += incr; 
                   
            theSys.particles[i].pos[k] = pos(k);
            theSys.particles[i].vel[k] = incr/dt;
        }
    }
    anGen->step(dt); //advance active noise in time
    theSys.apply_pbc();
    theSys.time++;
}

std::vector<arma::vec> ActiveSolver::get_active_noise_forces(System &theSys, Generator &gen)
{
    std::vector<arma::vec> active_forces(theSys.N);

    //Initialize to zero
    for (int i=0; i<theSys.N; i++) {
        arma::vec v(theSys.dim,arma::fill::zeros);
        active_forces[i] = v;
    } 

    //Compute real-space noise
    arma::field<arma::cx_vec> xi = gen.get_xi_r(1); //do_fft=1

    if (theSys.dim==3) {
        //Assign active force to each particle based on location in noise grid
        for(int i=0; i<theSys.N; i++) {
            arma::vec pos = theSys.particles[i].get_pos();
            int xind = std::round((pos(0)+0.5*theSys.edges[0])/gen.dx);
            int yind = std::round((pos(1)+0.5*theSys.edges[1])/gen.dx);
            int zind = std::round((pos(2)+0.5*theSys.edges[2])/gen.dx);
            for(int j=0; j<3; j++) active_forces[i](j) = va*xi(xind,yind,zind)(j).real(); //va is the active velocity
        }
    }
    else if (theSys.dim==2) {
        //Assign active force to each particle based on location in noise grid
        for(int i=0; i<theSys.N; i++) {
            arma::vec pos = theSys.particles[i].get_pos();
            int xind = std::round((pos(0)+0.5*theSys.edges[0])/gen.dx);
            int yind = std::round((pos(1)+0.5*theSys.edges[1])/gen.dx);
            //std::cout << i << "x comp: " << va*xi(xind,yind)(0).real() << std::endl;
            //std::cout << i << "y comp: " << va*xi(xind,yind)(1).real() << std::endl;
            for(int j=0; j<2; j++) active_forces[i](j) = va*xi(xind,yind)(j).real(); //va is the active velocity
        }
    }
    else if (theSys.dim==1) {
        //Assign active force to each particle based on location in noise grid
        for(int i=0; i<theSys.N; i++) {
            arma::vec pos = theSys.particles[i].get_pos();
            int xind = std::round((pos(0)+0.5*theSys.edges[0])/gen.dx);
            active_forces[i](0) = va*xi(xind)(0).real(); //va is the active velocity
        }
    }

    return active_forces;
}