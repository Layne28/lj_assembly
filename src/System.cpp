#include "System.hpp"

System::System(ParamDict &theParams, gsl_rng *&the_rg) {

    time = 0;

    //Set RNG
    rg = the_rg;

    //First assign default values to parameters
    kT = 1.0;
    rho = 4.0;
    phi = 0.4;
    dim = 3;
    dt = 0.005;

    obs = nullptr;

    particle_protocol = "zeros";

    sigma = 1.0;
    epsilon = 1.0;
    rcut = 2.5;

    //Then assign from ParamDict if there
    do_paramdict_assign(theParams);

    //Initialize Particles
    do_particle_init();

    //Initialize periodic image indices
    for (int i=0; i<N; i++) {
        std::vector<int> index(dim, 0);
        image.push_back(index);
    }
    this->apply_pbc();

    //Update old position to equal position after pbc
    for(int i=0; i<N; i++){
        particles[i].old_pos = particles[i].pos;
    }

    //Initialize neighbor grid
    if(do_neighbor_grid==1){
        std::array<double,2> min = {-0.5*edges[0], -0.5*edges[1]};
        std::array<double,2> max = {0.5*edges[0], 0.5*edges[1]};
        std::array<bool,2> p = {true, true};
        grid = new NeighborGrid<Particle, 2>(min, max, p, rcut);
        update_neighborgrid();
    }

    //Initialize cell list
    if(do_cell_list==1){
        if(dim!=2) {
            std::cout << "Error: cell lists not yet supported in d!=2." << std::endl;
            exit(-1);
        }
        ncell_x = int(edges[0]/rcut);
        ncell_y = int(edges[1]/rcut);
        cellsize_x = edges[0]/ncell_x;
        cellsize_y = edges[1]/ncell_y;
        if (fabs(ncell_x*cellsize_x-edges[0]>1e-3) || fabs(ncell_y*cellsize_y-edges[1]>1e-3)){
            std::cout << "Error: cell list dimensions computed incorrectly." << std::endl;
            exit(-1);
        }
        head = new int[ncell_x*ncell_y];
        list = new int[N];
        cellndx = new int[N];
        cellneigh = new int *[ncell_x*ncell_y];
        for(int i=0; i<ncell_x*ncell_y; i++) {
            cellneigh[i] = new int[9]; //at most 8 neighbor cells in 2D
        }
        fill_cellneigh();
    }
}

System::~System() {

    //clean up cell list
    if(do_cell_list==1){
        delete[] head;
        delete[] list;
        delete[] cellndx;
        for(int i=0; i<ncell_x*ncell_y; i++) {
            delete[] cellneigh[i];
        }
        delete[] cellneigh;
    }
}

void System::do_paramdict_assign(ParamDict &theParams) {

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

    if(theParams.is_key("kT")) kT = std::stod(theParams.get_value("kT"));
    if(theParams.is_key("phi")) phi = std::stod(theParams.get_value("phi"));
    if(theParams.is_key("dt")) dt = std::stod(theParams.get_value("dt"));
    if(theParams.is_key("particle_protocol")) particle_protocol = theParams.get_value("particle_protocol");
    if(theParams.is_key("epsilon")) epsilon = std::stod(theParams.get_value("epsilon"));
    if(theParams.is_key("sigma")) sigma = std::stod(theParams.get_value("sigma"));
    if(theParams.is_key("rcut")) rcut = std::stod(theParams.get_value("rcut"));
    if(theParams.is_key("do_cell_list")) do_cell_list = std::stoi(theParams.get_value("do_cell_list"));
    if(theParams.is_key("do_neighbor_grid")) do_neighbor_grid = std::stoi(theParams.get_value("do_neighbor_grid"));

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
    if (phi==0.0) N = 1;
    std::cout << "No. of particles: " << N << std::endl;
}

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

void System::do_particle_init() {

    if (particle_protocol=="zeros") {
        for (int i=0; i<N; i++) {
            Particle p(dim);
            particles.push_back(p);
        }
    }
    else if (particle_protocol=="random") {
        std::cout << "Doing random initialization..." << std::endl;
        //Assign particles random positions
        for (int i=0; i<N; i++) {
            int found = 0;
            while (found==0) {

                //Initialize a new particle at a random position
                Particle p(dim);
                for (int j=0; j<dim; j++) {
                    p.pos[j] = edges[j]*(gsl_rng_uniform(rg)-0.5);
                }
                
                found = 1;
                //Check for overlaps with existing particles
                for (int j=0; j<i; j++) {
                    if (get_dist(p, particles[j])<sigma*std::pow(2,1.0/6.0)) {
                        found = 0;
                    }
                }
                
                //If no overlaps found, we're done
                //Add particle to system
                if (found == 1) particles.push_back(p);
            }
        }
        std::cout << "Completed random initialization." << std::endl;
    }
    else if (particle_protocol=="lattice") {
        int nx = ceil(sqrt(N*sqrt(3.0)/2.0));
        int ny = ceil(sqrt(N*2.0/sqrt(3.0)));

        if (nx*ny<N){
            std::cout << "Error: number of lattice sites (" << nx*ny << ") is less than N." << std::endl;
            exit(-1);
        }
        int count = 0;
        double a = sqrt(edges[0]*edges[1]/((sqrt(3.0)*nx*ny/2)));
        std::cout << "spacing: " << a << std::endl;
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                if (count < N){
                    Particle p(2);
                    p.pos[0] = i*a;//-edges[0]/2.0;
                    if (j%2==1) p.pos[0] += 0.5*a;
                    p.pos[1] = (sqrt(3.0)/2.0)*j*a;//-edges[1]/2.0;
                    p.pos[2] = 0.0;
                    p.old_pos = p.pos;
                    count++;
                    particles.push_back(p);
                }
            }
        }
        if (count!=N){
            std::cout << "Error: didn't put as many particles on the lattice as should be there." << std::endl;
            exit(-1);
        }
        //Check for overlaps
        for(int i=0; i<N-1; i++){
            for(int j=i+1; j<N; j++){
                if (get_dist(particles[i], particles[j])<sigma*std::pow(2,1.0/6.0)) {
                    std::cout << get_dist(particles[i], particles[j]) << std::endl;
                    std::cout << "Warning: particles " << i << " and " << j << " are only separated by a distance " << get_dist(particles[i], particles[j]) << std::endl;
                }
            }
        }
    }
    else {
        throw std::runtime_error("Error: Particle initialization protocol not supported.");
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

double System::get_energy() {
    if (do_cell_list==1) return get_energy_cell_list();
    else if (do_neighbor_grid==1) return get_energy_neighbor_grid();
    else{
        double energy = 0;
        for (int i=0; i<N-1; i++) {
            for (int j=i+1; j<N; j++) {
                energy += get_energy_between(particles[i], particles[j]);
            }
        }
        return energy;
    }
}

double System::get_energy_between(Particle &p1, Particle &p2) {

    double energy = 0;
    double dist = get_dist(p1, p2);
    if (dist<=rcut) {
        if (potential_type=="lj"){
            energy = get_lj_potential(dist, sigma, epsilon, rcut);
        }
        else if(potential_type=="wca"){
            energy = get_wca_potential(dist, sigma, epsilon);
        }
    }

    return energy;
}

double System::get_energy_cell_list() {

    double energy = 0;

    //update cell list
    create_cell_list();

    //Loop over cells (including self)
    for (int index1 = N-1; index1 >= 0; index1--) {
        int icell = cellndx[index1];
        int index2 = head[icell];
        for (int nc = 0; nc < cellneigh[icell][0]; nc++) {
            int jcell = cellneigh[icell][nc+1];
            index2 = head[jcell];
            while (index2 != -1) {
                if(index1>index2){
                    double de = get_energy_between(particles[index1], particles[index2]);
                    energy += de;
                }
                index2 = list[index2];
            }
        }
    }
    return energy;//*0.5; //correct for double-counting
}

double System::get_energy_neighbor_grid() {

    double energy = 0;

    update_neighborgrid();

    for(int i=0; i<N; i++){
        /*
        for (auto p2 = grid->get_neighbor_iterator(&particles[i]); !p2.isDone(); p2++){ 
            std::cout << *p2 << std::endl;
        */
        std::vector<Particle *> neighbors = grid->get_neighbors(&particles[i]);
        for(unsigned int j=0; j<neighbors.size(); j++){
            Particle *p2 = neighbors[j];
            double dist = get_dist(particles[i], *p2);
            if (dist<1e-10) std::cout << "Error: particles overlap!" << std::endl;
            if (dist<=rcut) {
                if (potential_type=="lj"){
                    energy += get_lj_potential(dist, sigma, epsilon, rcut); 
                }
                else if(potential_type=="wca"){
                    energy += get_wca_potential(dist, sigma, epsilon); 
                }
            }
        }
    }
    return energy/2.0; //correct for double-counting neighbors
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

double System::get_wca_potential(double r, double sig, double eps) {

    double rc = sig*pow(2,1.0/6.0);
    if (r<rc) return get_lj_potential(r, sig, eps, rc);
    else return 0;
}

//O(N^2) calculation of forces (w/o cell list)

std::vector<arma::vec> System::get_forces() {
    
    if(do_cell_list==1) return get_forces_cell_list();
    else{
        std::vector<arma::vec> forces(N);
        for(int i=0; i<N; i++) {
            arma::vec force(dim);
            if(do_neighbor_grid==1) force = get_force_neighbor_grid(particles[i]);
            else force = get_force(particles[i]);
            forces[i] = force;
        }
        return forces;
    }
}


arma::vec System::get_force(Particle &p1) {
    
    arma::vec force(dim, arma::fill::zeros);
    
    for (int j=0; j<N; j++) {
        if (particles[j].get_id()!=p1.get_id()) {
            double dist = get_dist(p1, particles[j]);
            if (dist<1e-15) throw std::runtime_error("ERROR: attempting to divide by zero in force calculation!");
            if (dist<=rcut) {
                arma::vec disp = get_disp_vec(p1, particles[j]);
                if (potential_type=="lj"){
                    force += get_lj_force(dist, disp, sigma, epsilon); 
                }
                else if(potential_type=="wca"){
                    force += get_wca_force(dist, disp, sigma, epsilon); 
                }
            }
        }
    }
    return force;
}

arma::vec System::get_force_from(Particle &p1, Particle &p2) {
    
    arma::vec force(dim, arma::fill::zeros);
    
    double dist = get_dist(p1, p2);
    if (dist<1e-15) throw std::runtime_error("ERROR: attempting to divide by zero in force calculation!");
    if (dist<=rcut) {
        arma::vec disp = get_disp_vec(p1, p2);
        if (potential_type=="lj"){
            force = get_lj_force(dist, disp, sigma, epsilon); 
        }
        else if(potential_type=="wca"){
            force = get_wca_force(dist, disp, sigma, epsilon); 
        }
    }
    return force;
}

//Force with neighbor grid
arma::vec System::get_force_neighbor_grid(Particle &p1) {

    arma::vec force(dim, arma::fill::zeros);
    /*
    for (auto p2 = grid->get_neighbor_iterator(&p1); !p2.isDone(); p2++){
    */
    std::vector<Particle *> neighbors = grid->get_neighbors(&p1);
    for(unsigned int j=0; j<neighbors.size(); j++){
        Particle *p2 = neighbors[j];  
        double dist = get_dist(p1, *p2);
        if (dist<=rcut) {
            arma::vec disp = get_disp_vec(p1, *p2);
            if (potential_type=="lj"){
                force += get_lj_force(dist, disp, sigma, epsilon); 
            }
            else if(potential_type=="wca"){
                force += get_wca_force(dist, disp, sigma, epsilon); 
            }
        }
    }
    return force;
}

//TODO: fill this in

std::vector<arma::vec> System::get_forces_cell_list() {

    //update cell list
    create_cell_list();

    //Initialize forces to zero
    std::vector<arma::vec> forces(N);
    for(int i=0; i<N; i++){
        arma::vec force(dim, arma::fill::zeros);
        forces[i] = force;
    }

    //Loop over cells (including self)
    for (int index1 = N-1; index1 >= 0; index1--) {
        int icell = cellndx[index1];
        int index2 = head[icell];
        for (int nc = 0; nc < cellneigh[icell][0]; nc++) {
            int jcell = cellneigh[icell][nc+1];
            index2 = head[jcell];
            while (index2 != -1) {
                if(index1>index2){
                    arma::vec fij = get_force_from(particles[index1], particles[index2]);
                    forces[index1] += fij;
                    forces[index2] -= fij; //Newton's 3rd law
                }
                index2 = list[index2];
            }
        }
    }

    return forces;
}


arma::vec System::get_lj_force(double r, arma::vec rvec, double sig, double eps) {

    double rat = sig/r;
    double r2 = rat*rat;
    double r6 = r2*r2*r2;
    double r12 = r6*r6;

    double r14 = r12*r2;
    double r8 = r6*r2;

    return 24*(eps/sig)*(2*r14-r8)*rvec;
}

arma::vec System::get_wca_force(double r, arma::vec rvec, double sig, double eps) {

    arma::vec force = arma::zeros(dim);

    double rc = sig*pow(2,1.0/6.0);
    if (r<rc) return get_lj_force(r, rvec, sig, eps);
    else return force;
}

double System::minimize_energy(double tol) {
    std::cout << "Starting energy: " << get_energy() << std::endl;
    double diff_norm = 1.0; //normalized energy difference bt current and previous steps
    int max_iter = 1000;
    double eps = 1e-3;//1e-4; //gradient descent step size

    for(int t=0; t<max_iter; t++) {
        double e_old = get_energy();
        
        //Get forces
        std::vector<arma::vec> potential_forces = get_forces();

        //Update positions
        for (int i=0; i<N; i++) {
            particles[i].old_pos = particles[i].pos;
            for (int k=0; k<dim; k++) {
                particles[i].pos[k] += potential_forces[i](k)*eps;
            }
        }
        apply_pbc();
        double e_new = get_energy();
        diff_norm = fabs((e_new-e_old)/e_old);
        if (diff_norm<tol){
            std::cout << "Minimized energy after " << (t+1) << " iterations." << std::endl;
            std::cout << "End energy: " << e_new << std::endl;
            return e_new;
        } 
    }
    std::cout << "Attempted energy minimization for max iterations without converging." << std::endl;
    std::cout << "End energy: " << get_energy() << std::endl;
    return get_energy();
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

void System::update_neighborgrid() {

    for(int i=0; i<N; i++) {
        grid->update_atom(&particles[i]);
    }
}

/*** Cell list methods ***/

// assign each particle to a cell
void System::create_cell_list() {
   int icell;
   for (int i = 0; i < ncell_x*ncell_y; i++) {
      head[i] = -1;
   }
   for (int i = 0; i < N; i++) {
      //Assume pbc in [-L_mu/2, L_mu/2] for mu=x,y
      double shiftx = particles[i].pos[0] + edges[0]/2.0;
      double shifty = particles[i].pos[1] + edges[1]/2.0;
      icell = int(shiftx/cellsize_x) + int(shifty/cellsize_y)*ncell_x;
      if(icell>ncell_x*ncell_y) std::cout << "WARNING: icell=" << icell << std::endl;
      cellndx[i] = icell;
      list[i] = head[icell];
      if(list[i]>=N){
         std::cout << "WARNING: list[i]=" << list[i] << std::endl;
         exit(-1);
      }
      head[icell] = i;
   }
}

// find neighboring cells of each cell
void System::fill_cellneigh() {
    int icell, jcell, nneigh;
    int jx, jy;
    for (int ix = 0; ix < ncell_x; ix++) {
        for (int iy = 0; iy < ncell_y; iy++) {
            icell = ix + iy*ncell_x;
            nneigh = 0;
            for (int i = -1; i < 2; i++) {
                jx = ix+i;
                //Enforce pbc
                if(jx<0) jx += ncell_x;
                if(jx>=ncell_x) jx -= ncell_x;
                for (int j = -1; j < 2; j++) {
                    jy = iy+j;
                    //Enforce pbc
                    if(jy<0) jy += ncell_y;
                    if(jy>=ncell_y) jy -= ncell_y;
                    jcell = jx + jy*ncell_x;
                    cellneigh[icell][nneigh+1] = jcell;
                    nneigh++;
                }
            }
            cellneigh[icell][0] = nneigh;
            if(nneigh!=9){
                std::cout << "Error: number of neighbors (" << nneigh << ") should be 9 including cell itself." << std::endl;
                exit(-1);
            }
      }
   }
}