import numpy as np
import os
import mdtraj
from analysis import Trajectory, aligned_rmsd
import multiprocessing as mp

class GromacsTrajectory(Trajectory):

    """ Trajectory object from Gromacs files
    Keyword arguments:
    trajectory_filename -- filename of the trajectory, e.g. in the .xtc file format
    coordinate_filename -- filename of the folded / initial structure, e.g. in the .gro file format
    topology_filename -- filename of topology, e.g. in a .top file format
    """
    def __init__(self, trajectory_filename, coordinate_filename, topology_filename):
        super().__init__()
        self._read_coordinates(coordinate_filename, trajectory_filename)
        self._read_native_contacts(topology_filename)
        self._initial_conformation = mdtraj.load(coordinate_filename).xyz[0]
    
    def _read_coordinates(self, coordinate_filename, trajectory_filename):
        try:
            t = mdtraj.load(trajectory_filename, top=coordinate_filename)
            self._time_steps = t.time
            self._chains = t.xyz
        except MemoryError as e:
            trajectory_size = int(str(e).split('(')[1].split(')')[0].split(', ')[0])
            print(f'Warning: Using large trajectory ({trajectory_size} frames) mode with lower coordinate precision')
            t = mdtraj.load(coordinate_filename)
            self._chains = np.empty((trajectory_size, *t.xyz.shape[1:]), dtype=np.half)
            self._time_steps = np.empty((trajectory_size))
            with mdtraj.formats.XTCTrajectoryFile(trajectory_filename, mode='r') as file:
                step_size = 10000
                for i in range(trajectory_size // step_size + 1):
                    print(f'Reading frames {i*step_size} to {(i+1)*step_size-1} out of {trajectory_size}', end='\r')
                    xyz, time, step, box = file.read(step_size)
                    self._time_steps[i*step_size:(i+1)*step_size] = time
                    self._chains[i*step_size:(i+1)*step_size] = xyz.astype(np.half)
        self._residue_names = [i.residue.name for i in t.topology.atoms]
        self._box_dimensions = t.unitcell_lengths[0]
        residue_indices = np.array([i.residue.resSeq for i in t.topology.atoms])
        self._chain_indices = np.where(residue_indices == 1)[0]

    def _read_native_contacts(self, topology_filename):
        file = open(topology_filename, 'r')
        linestack = []
        status = 0
        non_bonded = []
        non_bonded_distances = []
        molecules = []
        atom_indices = {}
        atom_types = {}
        pairs = {}
        pairs_distances = {}
        name = ''
        inter_pairs = []
        inter_pairs_distances = []
        atom_types_all = []
        intra_pairs = []
        intra_pairs_distances = []
        count = 0
        for line in file:
            if '#include "' in line:
                subfilename = line.split('"')[1]
                filepath_list = topology_filename.split(os.sep)
                filepath_list[-1] = subfilename
                subfile = open(os.sep.join(filepath_list), 'r')
                linestack.extend(subfile.readlines())
                subfile.close()
            else:
                linestack.append(line)
            while len(linestack) > 0:
                data = linestack.pop(0)
                if data[0] == ';':
                    continue
                if '[ nonbond_params ]' in data:
                    status = 1
                elif '[ atoms ]' in data:
                    status = 2
                elif '[ pairs ]' in data:
                    status = 3
                elif '[ molecules ]' in data:
                    status = 4
                elif '[ moleculetype ]' in data:
                    status = 5
                elif len(data) < 2:
                    status = 0
                elif status == 1:
                    data_split = [x for x in data.split() if x]
                    non_bonded.append(data_split[:2])
                    non_bonded_distances.append(float(data_split[6].split(',')[0]))
                elif status == 5:
                    name = data.split()[0]
                    molecules.append(name)
                    pairs[name] = []
                    pairs_distances[name] = []
                    atom_indices[name] = []
                    atom_types[name] = []
                elif status == 2:
                    data_split = [x for x in data.split() if x]
                    atom_indices[name].append(int(data_split[0]))
                    atom_types[name].append(data_split[1])
                elif status == 3:
                    data_split = [x for x in data.split() if x]
                    pairs[name].append(list(map(int, data_split[:2])))
                    pairs_distances[name].append(float(data_split[6]))
                elif status == 4:
                    data_split = [x for x in data.split() if x]
                    m_name = data_split[0]
                    number = int(data_split[1])
                    for m in range(number):
                        for k in pairs[m_name]:
                            intra_pairs.append([atom_indices[m_name].index(k[0]) + count, atom_indices[m_name].index(k[1]) + count])
                        intra_pairs_distances.extend(pairs_distances[m_name])
                        atom_types_all.extend(atom_types[m_name])
                        count += len(atom_indices[m_name])
        atom_types_all = np.array(atom_types_all)
        for i, nb in enumerate(non_bonded):
            for i1 in np.where(atom_types_all == nb[0])[0]:
                for i2 in np.where(atom_types_all == nb[1])[0]:
                    inter_pairs.append([i1, i2])
                    inter_pairs_distances.append(non_bonded_distances[i])
        self._contacts_between_chains = np.array(inter_pairs)
        self._contacts_between_chains_distances = np.array(inter_pairs_distances)
        self._contacts_inside_chains = np.array(intra_pairs)
        self._contacts_inside_chains_distances = np.array(intra_pairs_distances)
        file.close()

    def write_coordinates(self, filename, chainId = None, title='Coordinate file'):
        if not self.is_valid:
            raise Error('Trajectory is invalid and can not be written')
        file = open(filename, 'w')
        bead_count = 0
        start_index, stop_index = self._get_chain_start_stop_indices(chainId)
        bead_number = stop_index - start_index
        for f in range(self._chains.shape[0]):
            file.write('{} t= {:.4f} step= {}\n  {}\n'.format(title, i*dt, i, bead_number))
            for i in range(start_index, stop_index):
                if i in self._chain_indices:
                    bead_count = 1
                file.write('{:>5}{:<5}{:>5}{:>5}{:>8}{:>8}{:>8}\n'.format(bead_count, self._residue_names[i], 'CA', bead_count, *['{:.3f}'.format(e) for e in self._chains[f, i]]))
                bead_count += 1
            file.write('{:>12} {:>12} {:>12}\n'.format(*['{:.7f}'.format(c) for c in self._box_dimensions]))
            bead_count = 1
        file.close()

class RMSD2D:

    @staticmethod
    def __init(shared_):
        global shared
        shared = shared_

    @staticmethod
    def __calculate_rmsd(indices):
        arr = np.frombuffer(shared[0], shared[1]).reshape(shared[2])
        for i, j in indices:
            arr[i, j] = arr[j, i] = aligned_rmsd(shared[3][i].xyz[0], shared[3][j].xyz[0])

    @staticmethod
    def calculate(trajectory_filename, coordinate_filename, step_size=1, processes=10):
        trajectory = mdtraj.load(trajectory_filename, top=coordinate_filename)[::step_size]
        size = trajectory.xyz.shape[0]
        shape = (size, size)
        data_type = np.dtype(float)
        shared_arr = mp.RawArray(np.ctypeslib.as_ctypes_type(data_type), size**2)
        with mp.Pool(processes, initializer=RMSD2D.__init, initargs=((shared_arr, data_type, shape, trajectory),)) as pool:
            indices = [(i, j) for i in range(size) for j in range(i + 1, size)]
            n = len(indices) // processes
            slices = [[n*i, n*(i+1)] for i in range(processes)]
            slices[-1][-1] = len(indices)
            pool.map(RMSD2D.__calculate_rmsd, [indices[s[0]:s[1]] for s in slices])
        arr = np.frombuffer(shared_arr, data_type).reshape(shape)
        return arr