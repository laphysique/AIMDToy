from pyscf import gto, scf, dft, grad

mp = 1836 

class AIMD(object):
    
    def __init__(self, atoms, m, x, v, dt, method='scf'):
        self.atoms = atoms
        self.inds = [x for x in range(len(atoms)) ]
        self.m = m
        self.x = x
        self.v = v
        self.dt = dt
        self.energy = 0
        self.a = np.zeros((len(atoms), 3))
       
        
        if method == 'scf':
            # Hatree-Fock
            self.scanner = gto.M().set(unit='Bohr', verbose=False).apply(scf.RHF).apply(grad.RHF)
        elif method == 'dft': 
            # DFT
            self.scanner = gto.M().set(unit='Bohr', verbose=False).apply(dft.RKS).apply(grad.RKS)
        else:
            raise ValueError('Wrong quantum chemistry method')
            
    def map_struct(self, i):
        return (self.atoms[i], self.x[i])

    def v_verlet(self): 
        self.x += self.v * self.dt + 0.5 * self.a * self.dt * self.dt
        self.v += 0.5 * self.a * self.dt  
    
        s = self.scanner.as_scanner()
        atom = map(self.map_struct, self.inds)
        self.energy, gradient = s(atom)

        self.a = np.multiply(-gradient.T, 1/self.m).T / mp
    
        self.v += 0.5 * self.a * self.dt

