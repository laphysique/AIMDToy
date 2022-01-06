from aimd import AIMD
import numpy as np
import time

io_unit_length = "Angstrom"
io_unit_time = "fs"

bohr = 0.529177210903  # Angstrom
dt_atom = 2.4188843265857e-2 # in units of fs


if __name__ == '__main__':
    
    """
    atoms = ['O', 'H', 'H']
    m_atoms = np.array([16, 1, 1])
    x_atoms = np.array([[0.000, 0.000, 0.000],
                        [0.000, 0.000, 1.000],
                        [0.943, 0.000, -0.333]]) / bohr
    """
 
    atoms = ['O', 'H', 'H', 'O', 'H', 'H']
    m_atoms = np.array([16, 1, 1, 16, 1, 1])
    x_atoms = np.array([[ 0.001,   0.000,   2.976],
                        [ 0.932,  -0.000,   3.306],
                        [ 0.000,  -0.000,   1.987],
                        [-0.001,  -0.000,  -0.001],
                        [-0.466,   0.808,  -0.329],
                        [-0.466,  -0.808,  -0.329]]) / bohr
    

    inds = np.zeros(len(atoms))
    v_atoms = np.zeros((len(atoms),3)) 
    a_atoms = np.zeros((len(atoms),3)) 

    e_tot = 0

    aimd = AIMD(atoms=atoms, m=m_atoms, x=x_atoms, v=v_atoms, dt=1/dt_atom)


    f = open('h2o_test.pdb', 'w')

    t0 = time.time()
    i = 0
    for t in range(1000):
        f.write('MODEL\n')

        for j in range(len(aimd.atoms)):

            row = ["ATOM ",
                   str(i).rjust(5),
                   aimd.atoms[j].rjust(2) + "  ",  # atom name
                   "   ",  # residue name
                   "A".ljust(3),  # chain ID.. will need to extend this one day.
                   "1".ljust(5),
                   f"{aimd.x[j, 0] * bohr:.3f}".rjust(7),
                   f"{aimd.x[j, 1] * bohr:.3f}".rjust(7),
                   f"{aimd.x[j, 2] * bohr:.3f}".rjust(7),
                   "1.00".rjust(5),
                   "\n"]    
            f.write(" ".join(row))

            i += 1
  
        f.write('ENDMDL\n')

        aimd.v_verlet()

    f.close()

    print('time spent: ', time.time()-t0)

