#Create input file given parameters

import argparse
import os

def main():

    parser = argparse.ArgumentParser(description='Write active noise input file with given input parameters.')
    parser.add_argument('indir',
                        help='Directory where input file will be written.')
    #System args
    parser.add_argument('--phi', default='0.4')
    parser.add_argument('--dim', default='2')
    parser.add_argument('--Lx', default='50.0')
    parser.add_argument('--Ly', default='50.0')
    parser.add_argument('--is_p_x', default='1')
    parser.add_argument('--is_p_y', default='1')
    parser.add_argument('--particle_protocol', default='random')
    parser.add_argument('--sigma', default='1.0')
    parser.add_argument('--epsilon', default='1.0')
    parser.add_argument('--kT', default='0.0')
    parser.add_argument('--do_cell_list', default='1')
    parser.add_argument('--do_neighbor_grid', default='0')

    #Solver args
    parser.add_argument('--dt', default='0.0002')
    parser.add_argument('--gamma', default='1.0')
    parser.add_argument('--va', default='0.1')

    #Observer args
    parser.add_argument('--output_dir', default='lj_assembly')
    parser.add_argument('--particles_freq', default='200')
    parser.add_argument('--thermo_freq', default='200')

    #LabBench args
    parser.add_argument('--experiment', default='standard')
    parser.add_argument('--equil_steps', default='50000')
    parser.add_argument('--production_steps', default='5000000')
    parser.add_argument('--info_freq', default='8000')

    #Active Noise Generator args
    parser.add_argument('--nx', default='100')
    parser.add_argument('--ny', default='100')
    parser.add_argument('--tau', default='1.0')
    parser.add_argument('--Lambda', default='1.0')

    args = parser.parse_args()
    print(args.__dict__)

    try:
        os.makedirs(args.indir)
        print('Made directory.')
    except OSError as e:
        print('Directory exists')

    print('Writing input file...')
    with open(args.indir + '/lj_assembly_kT=%.01f_phi=%.02f_va=%.02f_tau=%.02f_lambda=%.02f_Lx=%.01f_Ly=%.01f_nx=%d_ny=%d.in' % (float(args.kT), float(args.phi), float(args.va),float(args.tau),float(args.Lambda),float(args.Lx),float(args.Ly), int(args.nx), int(args.ny)), 'w') as f:

        f.write('#System\n')
        f.write('phi = %s\n' % args.phi)
        f.write('dim = %s\n' % args.dim)
        f.write('Lx = %s\n' % args.Lx)
        f.write('Ly = %s\n' % args.Ly)
        f.write('is_p_x = %s\n' % args.is_p_x)
        f.write('is_p_y = %s\n' % args.is_p_y)
        f.write('particle_protocol = %s\n' % args.particle_protocol)
        f.write('sigma = %s\n' % args.sigma)
        f.write('epsilon = %s\n' % args.epsilon)
        f.write('kT = %s\n' % args.kT)
        f.write('do_cell_list = %s\n' % args.do_cell_list)
        f.write('do_neighbor_grid = %s\n' % args.do_neighbor_grid)
        f.write('\n')

        f.write('#Solver\n')
        f.write('dt = %s\n' % args.dt)
        f.write('gamma = %s\n' % args.gamma)
        f.write('va = %s\n' % args.va)
        f.write('\n')

        f.write('#Observer\n')
        f.write('output_dir = %s\n' % args.output_dir)
        f.write('particles_freq = %s\n' % args.particles_freq)
        f.write('thermo_freq = %s\n' % args.thermo_freq)
        f.write('\n')

        f.write('#LabBench\n')
        f.write('experiment = %s\n' % args.experiment)
        f.write('equil_steps = %s\n' % args.equil_steps)
        f.write('production_steps = %s\n' % args.production_steps)
        f.write('info_freq = %s\n' % args.info_freq)
        f.write('\n')

        f.write('#Active Noise Generator\n')
        f.write('nx = %s\n' % args.nx)
        f.write('ny = %s\n' % args.ny)
        f.write('tau = %s\n' % args.tau)
        f.write('lambda = %s\n' % args.Lambda)

if __name__=='__main__':
    main()
