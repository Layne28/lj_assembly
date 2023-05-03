#include "System.hpp"

System::System(ParamDict &theParams) {

    time = 0;

    //First assign default values to parameters
    kT = 1.0;
    rho = 4.0;
    phi = 0.4;
    dim = 3;
    dt = 0.005;

    obs = nullptr;

    //These are temporary, won't need to be accessed again outside of constructor
    std::string particle_protocol = "zeros";

    sigma = 1.0;
    epsilon = 1.0;
    rcut = 2.5;

    //Get number of dimensions from ParamDict
    if(theParams.is_key("dim")) dim = std::stoi(theParams.get_value("dim"));
    //Make "edges" and "is_periodic" vectors the correct size
    for(int i=0; i<dim; i++) {
        //Get box edges
        edges.push_back(10.0); //default value
        if(i==0 && theParams.is_key("Lx")) edges[i] = std::stod(theParams.get_value("Lx")); //read in value
        if(i==1 && theParams.is_key("Ly")) edges[i] = std::stod(theParams.get_value("Ly")); //read in value
        if(i==2 && theParams.is_key("Lz")) edges[i] = std::stod(theParams.get_value("Lz")); //read in value

        //Get whether each dimension is periodic
        is_periodic.push_back(1); //default value
        if(i==0 && theParams.is_key("is_p_x")) is_periodic[i] = std::stoi(theParams.get_value("is_p_x")); //read in value
        if(i==1 && theParams.is_key("is_p_y")) is_periodic[i] = std::stoi(theParams.get_value("is_p_y")); //read in value
        if(i==2 && theParams.is_key("is_p_z")) is_periodic[i] = std::stoi(theParams.get_value("is_p_z")); //read in value
    }

    //Then assign from ParamDict if there
    if(theParams.is_key("kT")) kT = std::stod(theParams.get_value("kT"));
    if(theParams.is_key("phi")) phi = std::stod(theParams.get_value("phi"));
    if(theParams.is_key("dt")) dt = std::stod(theParams.get_value("dt"));
    if(theParams.is_key("particle_protocol")) particle_protocol = theParams.get_value("particle_protocol");
    if(theParams.is_key("epsilon")) epsilon = std::stod(theParams.get_value("epsilon"));
    if(theParams.is_key("sigma")) sigma = std::stod(theParams.get_value("sigma"));
    if(theParams.is_key("rcut")) rcut = std::stod(theParams.get_value("rcut"));

    //Compute parameters that are inferred from other parameters
    double volume = 1.0;
    for(int i=0; i<dim; i++){
        volume *= edges[i];
    }
    double particle_volume;
    if (dim==1) {
        particle_volume = sigma;
    }
    else if (dim==2) {
        particle_volume = M_PI*(sigma/2)*(sigma/2);
    }
    else if (dim==3) {
        particle_volume = (4.0/3.0)*M_PI*(sigma/2)*(sigma/2)*(sigma/2);
    }
    else {
        std::cout << "Error: " << dim << "-dimensional simulation not supported." << std::endl;
        exit(-1);
    }
    //Compute N from volume and packing fraction
    N = std::round(phi*volume/particle_volume); 
    std::cout << "No. of particles: " << N << std::endl;

    //Initialize Particles
    if (particle_protocol=="zeros") {
        for (int i=0; i<N; i++) {
            Particle p(dim);
            particles.push_back(p);
        }
    }
    else if (particle_protocol=="random") {
        std::cout << "to be implemented" << std::endl;
    }
    else {
        throw std::runtime_error("Error: Particle initialization protocol not supported.");
    }

    this->apply_pbc();
    this->zero_com();

    //Initialize periodic image indices
    for (int i=0; i<N; i++) {
        std::vector<int> index(dim, 0);
        image.push_back(index);
    }
    
}

System::~System() {}

void System::apply_pbc() {

    for (int i=0; i<N; i++) {
        for (int j=0; j<dim; j++) {
            if (is_periodic[j]) {
                if (particles[i].pos[j] < -0.5*edges[j]) {
                    particles[i].pos[j] += edges[j];
                    image[i][j] -= 1;
                }
                if (particles[i].pos[j] >= 0.5*edges[j]) {
                    particles[i].pos[j] -= edges[j];
                    image[i][j] += 1;
                }
            }
        }
    }
}

void System::zero_com() {

    arma::vec com(dim, arma::fill::zeros);
    for (int i=0; i<N; i++) {
        com += particles[i].pos;
    }
    for (int i=0; i<N; i++) {
        particles[i].pos -= com/N;
    }
}

arma::vec System::get_com() {

    arma::vec com(dim, arma::fill::zeros);
    for (int i=0; i<N; i++) {
        com += particles[i].pos;
    }
    return com;
}

void System::set_obs(Observer &anObs) {

    obs = &anObs;
}

//TODO: change this to work for arbitrary potential
double System::get_energy()
{
    double energy = 0;
    for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
            double dist = get_dist(particles[i], particles[j]);
            if (dist<=rcut) {
                energy += get_lj_potential(dist, sigma, epsilon, rcut);
            }
        }
    }
    return energy;
}

double System::get_lj_potential(double r, double sig, double eps, double rc) {

    double rat = sig/r;
    double r2 = rat*rat;
    double r6 = r2*r2*r2;
    double r12 = r6*r6;

    double rat_cut = sig/rc;
    double r2_cut = rat_cut*rat_cut;
    double r6_cut = r2_cut*r2_cut*r2_cut;
    double r12_cut = r6_cut*r6_cut;

    return 4*eps*(r12-r6) - 4*eps*(r12_cut-r6_cut);
}

//TODO: modify this so it calls a potential function (so we can switch out different potentials easily)
arma::vec System::get_force (Particle &p1) {

    arma::vec force(dim, arma::fill::zeros);
    
    for (int j=0; j<N; j++) {
        if (particles[j].get_id()!=p1.get_id()) {
            double dist = get_dist(p1, particles[j]);
            if (dist<1e-15) throw std::runtime_error("ERROR: attempting to divide by zero in force calculation!");
            if (dist<=rcut) {
                arma::vec disp = get_disp_vec(p1, particles[j]);
                force += get_lj_force(dist, disp, sigma, epsilon); 
            }
        }
    }
    return force;
}

arma::vec System::get_lj_force(double r, arma::vec rvec, double sig, double eps) {

    arma::vec force = arma::zeros(dim);

    double rat = sig/r;
    double r2 = rat*rat;
    double r6 = r2*r2*r2;
    double r12 = r6*r6;

    double r14 = r12*r2;
    double r8 = r6*r2;

    return 24*(eps/sig)*(2*r14-r8)*rvec;
}

arma::vec System::get_disp_vec(Particle &p1, Particle &p2) {

    arma::vec disp = arma::zeros(dim);
    for (int j=0; j<dim; j++) disp[j] = p1.pos[j]-p2.pos[j];
    for (int j=0; j<dim; j++) {
        if (is_periodic[j]) {
            if (disp[j]<(-0.5*edges[j])) disp[j] += edges[j];
            if (disp[j]>=(0.5*edges[j])) disp[j] -= edges[j];
        }
    }

    return disp;
}

double System::get_dist(Particle &p1, Particle &p2) {

    arma::vec disp = get_disp_vec(p1, p2);
    double len = 0;
    for (int i=0; i<dim; i++) len += disp[i]*disp[i];
    len = sqrt(len);

    return len;
}

Observer System::get_obs() {

    if (obs) {
        return *obs;
    }
    else {
        throw("Error: no observer set!");
    }
}