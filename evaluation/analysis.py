from re import S
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist
from scipy.special import expit
from scipy import optimize
import networkx as nx
import itertools
from process_helper import parmap
import matplotlib.pyplot as plt
from typing import List, Tuple

class Trajectory:

    def __init__(self) -> None:
        self._chains = np.array([])
        self._initial_conformation = np.array([])
        self._chain_indices = []
        self._residue_names = []
        self._time_steps = np.array([])
        self._box_dimensions = np.array([0, 0, 0])
        self._contacts_inside_chains = np.array([])
        self._contacts_between_chains = np.array([])
        self._contacts_inside_chains_distances = np.array([])
        self._contacts_between_chains_distances = np.array([])

    def __len__(self) -> int:
        return len(self._chains)

    def size(self, chainId: int=None) -> int:
        """ Get size of chain at index `chainId` """
        start_index, stop_index = self.get_chain_start_stop_indices(chainId)
        return stop_index - start_index

    def get_number_of_chains(self) -> int:
        """ Get the number of chains """
        return len(self._chain_indices)

    def get_trajectory(self, chainId: int=None) -> np.ndarray:
        """ Get the trajectory of the chain at index `chainId` """
        start_index, stop_index = self.get_chain_start_stop_indices(chainId)
        return self._chains[:,start_index:stop_index]

    def get_initial_conformation(self, chainId: int=None) -> np.ndarray:
        start_index, stop_index = self.get_chain_start_stop_indices(chainId)
        return self._initial_conformation[start_index:stop_index]

    def get_time_steps(self) -> np.ndarray:
        """ Get the timestamps in picoseconds corresponding to the frames in the trajectory """
        return self._time_steps

    def get_residue_names(self, chainId: int=None) -> List[str]:
        """ Get the 3-letter names of all residues in chain with index `chainId` """
        start_index, stop_index = self.get_chain_start_stop_indices(chainId)
        return self._residue_names[start_index:stop_index]

    def get_box_dimensions(self) -> np.ndarray:
        """ Get box dimensions in nanometers """
        return self._box_dimensions

    def get_native_pairs_inside_chains(self) -> np.ndarray:
        """ List of amino acid index pairs which are defined as native pairs within all chains """
        return self._contacts_inside_chains

    def get_native_pairs_between_chains(self) -> np.ndarray:
        """ List of amino acid index pairs which are defined as native pairs between chains """
        return self._contacts_between_chains

    def get_chain_start_stop_indices(self, chainId) -> Tuple[int, int]:
        """ Get the amino acid start and end indices of the chain at index `chainId` within the entire structure """
        if chainId == None:
            return 0, self._chains.shape[1]
        start_index = self._chain_indices[chainId]
        stop_index = self._chains.shape[1] if (chainId == len(self._chain_indices) - 1) else self._chain_indices[chainId + 1]
        return start_index, stop_index

    def __getitem__(self, pos):
        if type(pos) is int:
            start_index, stop_index = self.get_chain_start_stop_indices(pos)
            return self._chains[:, start_index:stop_index]
        elif type(pos) is tuple and len(pos) == 2:
            start_index, stop_index = self.get_chain_start_stop_indices(pos[0])
            return self._chains[pos[1], start_index:stop_index]
        elif type(pos) is tuple and len(pos) == 3:
            start_index, stop_index = self.get_chain_start_stop_indices(pos[0])
            return self._chains[pos[1], start_index:stop_index][pos[2]]

    def is_valid(self):
        """ Check if the trajectory object is valid """
        if len(self._chains) == 0 or len(self._chain_indices) == 0 or len(self._residue_names) == 0:
            return False
        if len(self._time_steps) == 0 or len(self._box_dimensions) != 3:
            return False
        if self._chains is not np.ndarray or self.chains.shape[1] != len(self._residue_names):
            return False
        if self._chains.shape[0] != len(self._time_steps):
            return False
        if len(self._chains.shape) != 3 or self._chains.shape[2] != 3:
            return False
        return True

    @staticmethod
    def __frame_to_slice(frame):
        if type(frame) is slice:
            return frame
        if type(frame) is tuple and len(frame) == 2:
            return slice(*frame, 1)
        if type(frame) is tuple and len(frame) == 3:
            return slice(*frame)
        if type(frame) is int:
            return slice(frame, frame + 1, 1)
        return slice(None)

    def _get_single_fraction_of_native_contact(self, coordinates: np.ndarray, native_pairs: np.ndarray, cut_off: float, sigmoid_factor: float, pbc: bool) -> int:
        """ Calculates the fraction of native contacts for a single trajectory frame for a given list of native pairs """
        a1 = coordinates[native_pairs[:,0]]
        a2 = coordinates[native_pairs[:,1]]
        if pbc:
            delta = np.abs(a1 - a2)
            distances = np.linalg.norm(np.where(delta > 0.5 * self._box_dimensions, delta - self._box_dimensions, delta), axis=1)
        else:
            distances = np.linalg.norm(a1 - a2, axis=1)
        if sigmoid_factor <= 0:
            return np.mean(distances < cut_off)
        else:
            return np.mean(expit(sigmoid_factor*(cut_off-distances)))
    
    def _get_fraction_of_native_contacts(self, coordinates: np.ndarray, native_pairs: np.ndarray, cut_off: float, sigmoid_factor: float, pbc: bool) -> np.ndarray:
        """ Calculates the fraction of native contacts for multiple frames for a given list of native pairs in parallel """
        if coordinates.shape[0] == 1:
            return self._get_single_fraction_of_native_contact(coordinates[0], native_pairs, cut_off, sigmoid_factor, pbc)
        else:
            result = []
            result = np.array(parmap(
                lambda frame: self._get_single_fraction_of_native_contact(frame, native_pairs, cut_off, sigmoid_factor, pbc), coordinates))
            return np.array(result)


    """ Get the fraction of native contacts for all intramolecular native pairs

    Keyword arguments:
    frame -- If provided, the a subset of frames from the trajectory is used
        int: single frame
        Tuple[int, int]: range of frames from start-index to end-index
        Tuple[int, int, int]: range of frames from start-index to end-index with stepsize
    chain_id -- If provided, onle the chain at this index is included
    offset -- [nm] This value is added to the native pair distance to calculate the FNC cutoff
    sigmoid_factor -- If > 0, a sigmoidal cutoff is used instead of a hard cutoff
    pbc -- Consider periodic boundaries by calculating closest-image distances instead of real distances
    """
    def get_fraction_of_native_contacts_inside(self, frame=None, chain_id=None, offset: float=0.15, sigmoid_factor: float=20, pbc: bool=False) -> np.ndarray:
        frame_selection = Trajectory.__frame_to_slice(frame)
        chain_start, chain_stop = self.get_chain_start_stop_indices(chain_id)
        selector = (self._contacts_inside_chains >= chain_start).all(axis=1) & (self._contacts_inside_chains < chain_stop).all(axis=1)
        contacts = self._contacts_inside_chains[selector]
        cut_off = self._contacts_inside_chains_distances[selector] + offset
        return self._get_fraction_of_native_contacts(self._chains[frame_selection], contacts, cut_off, sigmoid_factor, pbc)

    """ Get the fraction of native contacts for all intermolecular native pairs

    Keyword arguments:
    frame -- If provided, the a subset of frames from the trajectory is used
        int: single frame
        Tuple[int, int]: range of frames from start-index to end-index
        Tuple[int, int, int]: range of frames from start-index to end-index with stepsize
    offset -- [nm] This value is added to the native pair distance to calculate the FNC cutoff
    sigmoid_factor -- If > 0, a sigmoidal cutoff is used instead of a hard cutoff
    pbc -- Consider periodic boundaries by calculating closest-image distances instead of real distances
    """
    def get_fraction_of_native_contacts_between(self, frame=None, offset: float=0.15, sigmoid_factor: float=20, pbc: bool=False) -> np.ndarray:
        frame_selection = Trajectory.__frame_to_slice(frame)
        digitized = np.digitize(self._contacts_between_chains, self._chain_indices) - 1
        contacts = self._contacts_between_chains[digitized[:,0] != digitized[:,1]]
        cut_off = self._contacts_between_chains_distances[digitized[:,0] != digitized[:,1]] + offset
        return self._get_fraction_of_native_contacts(self._chains[frame_selection], contacts, cut_off, sigmoid_factor, pbc)

    """ Calculate the fraction of native contacts for all possible chain pairs

    Keyword arguments:
    frame -- If provided, the a subset of frames from the trajectory is used
        int: single frame
        Tuple[int, int]: range of frames from start-index to end-index
        Tuple[int, int, int]: range of frames from start-index to end-index with stepsize
    offset -- [nm] This value is added to the native pair distance to calculate the FNC cutoff
    sigmoid_factor -- If > 0, a sigmoidal cutoff is used instead of a hard cutoff
    pbc -- Consider periodic boundaries by calculating closest-image distances instead of real distances
    """
    def get_fraction_of_native_contacts_chain_pairs(self, frame=None, offset: float=0.15, sigmoid_factor: float=20, pbc: bool=False) -> Tuple[np.ndarray, np.ndarray]:
        frame_selection = Trajectory.__frame_to_slice(frame)
        digitized = np.digitize(self._contacts_between_chains, self._chain_indices) - 1
        groups = np.unique(digitized, axis=0)
        result = []
        for group in groups:
            contacts = self._contacts_between_chains[(digitized[:,0] == group[0]) & (digitized[:,1] == group[1])]
            cut_off = self._contacts_inside_chains_distances[(digitized[:,0] == group[0]) & (digitized[:,1] == group[1])] + offset
            result.append(self._get_fraction_of_native_contacts(self._chains[frame_selection], contacts, cut_off, sigmoid_factor, pbc))
        return result, groups

    def _get_chain_slices(self, overwrite: dict=None):
        indices = [*self._chain_indices, len(self._residue_names)]
        slices = [slice(indices[i], indices[i+1], 1) for i in range(len(self._chain_indices))]
        if overwrite is not None and type(overwrite) == dict:
            for key, value in overwrite.items():
                slices[key] = slice(slices[key].start + value.start, slices[key].start + value.stop, value.step if value.step else 1)
        return slices

    """ Generate 2D distance matrices for distances between the geometrical centers of all chains

    Keyword arguments:
    pbc -- Consider periodic boundaries by calculating closest-image distances instead of real distances
    overwrite_chain_slices -- set the amino acid selection slice (value) for chain at index (key)
    """
    def get_chain_distance_map(self, pbc: bool=False, overwrite_chain_slices: dict=None) -> Tuple[np.ndarray, np.ndarray]:
        slices = self._get_chain_slices(overwrite=overwrite_chain_slices)
        number_of_chains = self.get_number_of_chains()
        init = np.zeros((number_of_chains, number_of_chains))
        traj = np.zeros((len(self._chains), number_of_chains, number_of_chains))
        for i in range(number_of_chains):
            for j in range(i + 1, number_of_chains):
                if pbc:
                    init_delta = np.abs(self._initial_conformation[slices[i]].mean(axis=0) - self._initial_conformation[slices[j]].mean(axis=0))
                    init[i,j] = init[j,i] = np.linalg.norm(np.where(init_delta > 0.5 * self._box_dimensions, init_delta - self._box_dimensions, init_delta))
                    traj_delta = np.abs(self._chains[:,slices[i]].mean(axis=1) - self._chains[:,slices[j]].mean(axis=1))
                    traj[:,i,j] = traj[:,j,i] = np.linalg.norm(np.where(traj_delta > 0.5 * self._box_dimensions, traj_delta - self._box_dimensions, traj_delta), axis=1)
                else:
                    init[i,j] = init[j,i] = np.linalg.norm(self._initial_conformation[slices[i]].mean(axis=0) - self._initial_conformation[slices[j]].mean(axis=0))
                    traj[:,i,j] = traj[:,j,i] = np.linalg.norm(self._chains[:,slices[i]].mean(axis=1) - self._chains[:,slices[j]].mean(axis=1), axis=-1)
        return init, traj

    """ Calculates pairwise chain distances over the entire trajectory

    Keyword arguments:
    rolling_window_size -- apply a rolling window averaging with the given size over the chain distances
    overwrite_chain_slices -- set the amino acid selection slice (value) for chain at index (key)
    processes -- Number of processes for parallelized distance calculations
    """
    def get_chain_distances(self, rolling_window_size: int=40, overwrite_chain_slices: dict=None, processes: int=None) -> np.ndarray:
        number_of_chains = self.get_number_of_chains()
        slices = self._get_chain_slices(overwrite=overwrite_chain_slices)
        indices = np.array([(i, j) for i in range(number_of_chains) for j in range(i + 1, number_of_chains)])
        number_of_frames = self._chains.shape[0]
        def calculate_distances(index):
            i, j = index
            distance = np.absolute(self._chains[:,slices[i]].mean(axis=1) - self._chains[:,slices[j]].mean(axis=1))
            distance = np.minimum(self._box_dimensions - distance, distance)
            distance = np.linalg.norm(distance, axis=-1)
            strides = distance.strides + (distance.strides[-1],)
            distance = np.lib.stride_tricks.as_strided(distance, shape=(number_of_frames - rolling_window_size + 1, rolling_window_size), strides=strides)
            return distance
        result = np.array(parmap(calculate_distances, indices, nprocs=processes))
        return result

    """ Calculate one graph representation for each frame in the trajectory
    In the graph, each chain is represented by one node. Two nodes are connected if their chains are bound.
    Whether two chains are bound is determined by their distance (in comparison to their native distance) and
    the variations in their distance.

    Keyword arguments:
    rolling_window_size -- apply a rolling window averaging with the given size over the chain distances.
        This reduces false positives of bound chains.
    std_cut_off -- cutoff value for the standard deviation of the chain distance between two bound chains
    distance_cutoff_factor -- a distance histogram is used to determine the native distance between chains.
        Increase/decrease this factor to allow for larger/smaller distance cutoffs.
    single_distance_cutoff -- If the histogram method fails (usually for small complexes), this factor is multiplied
        with the average chain distance to determine the distance cutoff
    fixed_distance_cutoff -- If specified, this value is used as a distance cutoff (instead of an automatic cutoff determination)
    overwrite_chain_slices -- set the amino acid selection slice (value) for chain at index (key)
    logging -- Print logging information
    show_distance_histogram -- Visualize the distance histogram used for the determination of the native chain distances
    processes -- Number of processes for parallelized chain contact calculations
    """
    def get_chain_graph(self, rolling_window_size: int=30, std_cut_off: float=0.3, distance_cutoff_factor: float=0.5, single_distance_cutoff: float=1.2,
            fixed_distance_cutoff: float=None, overwrite_chain_slices: dict=None, logging: bool=True, show_distance_histogram: bool=False, processes: int=None) -> List[nx.Graph]:
        number_of_chains = self.get_number_of_chains()
        original_slices = self._get_chain_slices()
        slices = self._get_chain_slices(overwrite=overwrite_chain_slices)
        # Get chains
        chains = {}
        for i in range(len(self._chain_indices)):
            seq = ''.join(self._residue_names[original_slices[i]])
            if seq in chains:
                chains[seq].append(i)
            else:
                chains[seq] = [i]
        keys = list(chains.keys())
        for i in range(len(keys)):
            chains[i] = chains[keys[i]]
            del chains[keys[i]]
        lookup = {}
        for key, indices in chains.items():
            for index in indices:
                lookup[index] = key
        # Get distances from native conformation
        distances = {}
        for name1 in chains.keys():
            for name2 in chains.keys():
                key = (name1, name2) if name1 < name2 else (name2, name1)
                if key not in distances:
                    distances[key] = []
                    for index1 in chains[name1]:
                        for index2 in chains[name2]:
                            if name1 != name2 or index1 > index2:
                                distances[key].append(np.linalg.norm(self._initial_conformation[slices[index1]].mean(axis=0) - self._initial_conformation[slices[index2]].mean(axis=0)))
                    if len(distances[key]) == 0:
                        del distances[key]
        for key in distances.keys():
            if len(distances[key]) <= len(self._chain_indices) * 2:
                short_distances = distances[key]
            else:
                short_distances = np.partition(distances[key], len(self._chain_indices) * 2)[:len(self._chain_indices) * 2]
            distance_hist = np.histogram(short_distances, bins=15)
            if show_distance_histogram:
                plt.figure()
                counts, bins = distance_hist
                plt.hist(bins[:-1], bins, weights=counts)
                plt.xlabel('Distance [nm]')
                plt.ylabel('Counts')
                plt.grid(alpha=0.2)
                plt.show()
            if fixed_distance_cutoff is not None:
                cutoff = fixed_distance_cutoff
            elif np.sum(distance_hist[0] != 0) > 1 and distance_hist[0].max() > 1 and distance_cutoff_factor < 1:
                distance_indices = np.where(distance_hist[0] > distance_hist[0].max() / 2)[0]
                cutoff = (distance_hist[1][distance_indices[0]+1] * distance_cutoff_factor + (1-distance_cutoff_factor) * distance_hist[1][distance_indices[1]])
            else:
                if logging:
                    print('Histogram method for the determination of native pair distances failed. Using average chain pair distances...')
                cutoff = np.array(distances[key]).mean() * single_distance_cutoff
            distances[key] = cutoff
        # Process trajectory
        indices = np.array([(i, j) for i in range(number_of_chains) for j in range(i + 1, number_of_chains)])
        number_of_frames = self._chains.shape[0]
        if logging:
            print('Calculating neighbors')
        def calculate_contacts(index):
            i, j = index
            distance = np.absolute(self._chains[:,slices[i]].mean(axis=1) - self._chains[:,slices[j]].mean(axis=1))
            distance = np.minimum(self._box_dimensions - distance, distance)
            distance = np.linalg.norm(distance, axis=-1)
            strides = distance.strides + (distance.strides[-1],)
            distance = np.lib.stride_tricks.as_strided(distance, shape=(number_of_frames - rolling_window_size + 1, rolling_window_size), strides=strides)
            return np.logical_and(distance.std(axis=-1) < std_cut_off, distance.mean(axis=-1) < distances[(lookup[i], lookup[j])])
        contacts = np.array(parmap(calculate_contacts, indices, nprocs=processes))
        if logging:
            print('Generating graphs')
        def generate_graph(frame):
            return nx.from_edgelist(indices[frame])
        graphs = parmap(generate_graph, contacts.T,  nprocs=processes)
        return graphs

    """ Calculate a clustered representation of the complex assembly process from a list of graph representations
    Only the number number and size of connected graph components (i.e. subassemblies) is used for clustering.
    This results in a smaller number of assembly states than generated by the `get_chain_formations_detailed`
    method. However, usually the smaller number of states is easier to interpret.

    Keyword arguments:
    graphs -- List of nx.Graph
    threshold -- All assembly states which appear less than this number of frames is excluded
    filtering -- Remove short appearances of assembly states via an moving average filtering with a window size of `filtering`
    logging -- Print logging information
    """
    def get_chain_formations(self, graphs, threshold: int=100, filtering: int=0, logging: bool=False):
        if logging:
            print('Get graph properties')
        components = []
        for graph in graphs:
            components.append(sorted([len(sg) for sg in nx.connected_components(graph)])[::-1])
        max_components = max([len(d) for d in components])
        for i in range(len(components)):
            while len(components[i]) < max_components:
                components[i].append(0)
        components = np.array(components)

        unique_labels = ['c{}'.format(i) for i in range(max_components)]

        if threshold > 0:
            if logging:
                print('Remove rare entries')
            unique, counts = np.unique(components, return_counts=True, axis=0)
            rare = unique[counts < threshold][::-1]
            unique = unique[::-1]
            if len(rare) > 0:
                if (rare[0] == unique[0]).all():
                    rare = rare[1:]
                rare_indices = np.nonzero((components[:, None] == rare).all(axis=-1).any(axis=-1))[0]
                for k, ri in enumerate(rare_indices):
                    nnind = np.nonzero(rare_indices[k:] - np.arange(ri, ri + len(rare_indices[k:]), 1))[0]
                    if len(nnind) == 0:
                        nnind = np.nonzero(rare_indices[:k+1] - np.arange(ri - k, ri+1, 1))[0]
                        if len(nnind) == 0:
                            continue
                        else:
                            index = ri - nnind[-1]
                    else:
                        index = ri + nnind[0]
                    components[ri] = components[index]

        unique = np.unique(components, axis=0)[::-1]
        img = np.equal(components[:,np.newaxis], unique).all(axis=2).T

        if filtering > 0:
            if logging:
                print('Filter result')
            filtered = []
            length = img.shape[1]
            for i in range(length):
                si_length = min(length, i+filtering) - max(0, i-filtering)
                si = img[:,max(0, i-filtering):min(length, i+filtering)].sum(axis=1)
                si[si < si_length * 0.5] = 0
                si[si > 0] = 1
                filtered.append(si)
            img = np.array(filtered).T
            zero_indices = (img.sum(axis=-1) > threshold)
            img = img[zero_indices]
            unique = np.array(unique)[zero_indices]

        if logging:
            print('Create labels')
        unique = unique.astype('<U8')
        zero_indices = ~np.all(unique == '0', axis=0)
        unique = unique[:, zero_indices]
        unique = [' '.join(u) for u in unique]

        return img, unique

    """ Calculate a clustered representation of the complex assembly process from a list of graph representations
    The number number and size of connected graph components (i.e. subassemblies) and the degree histogram of
    the graph is used for clustering. This results in more assembly states than generated by the `get_chain_formations`
    method.

    Keyword arguments:
    graphs -- List of nx.Graph
    threshold -- All assembly states which appear less than this number of frames is excluded
    filtering -- Remove short appearances of assembly states via an moving average filtering with a window size of `filtering`
    logging -- Print logging information
    """
    def get_chain_formations_detailed(self, graphs, threshold=100, filtering=0, logging=False):
        if logging:
            print('Get graph properties')
        degrees = []
        components = []
        for graph in graphs:
            degrees.append(nx.degree_histogram(graph) if len(graph.edges) > 0 else [])
            components.append(sorted([len(sg) for sg in nx.connected_components(graph)])[::-1])
        max_degree = max([len(d) for d in degrees])
        max_components = max([len(d) for d in components])
        for i in range(len(degrees)):
            while len(degrees[i]) < max_degree:
                degrees[i].append(0)
            degrees[i].extend(components[i])
            while len(degrees[i]) < max_degree + max_components:
                degrees[i].append(0)
        degrees = np.array(degrees)
        max_values = degrees.max(axis=0)

        unique_labels = ['d{}'.format(i) for i in range(max_degree)]
        unique_labels.append('|')
        unique_labels.extend(['c{}'.format(i) for i in range(max_components)])

        if threshold > 0:
            if logging:
                print('Remove rare entries')
            unique, counts = np.unique(degrees, return_counts=True, axis=0)
            rare = unique[counts < threshold]
            ind = np.lexsort([unique[:,i] for i in range(max_values[:max_degree].argsort()[-1] + 1)])
            unique = unique[ind][::-1]
            ind = np.lexsort([rare[:,i] for i in range(max_values[:max_degree].argsort()[-1] + 1)])
            rare = rare[ind][::-1]
            if len(rare) > 0:
                if (rare[0] == unique[0]).all():
                    rare = rare[1:]
                rare_indices = np.nonzero((degrees[:, None] == rare).all(axis=-1).any(axis=-1))[0]
                for k, ri in enumerate(rare_indices):
                    nnind = np.nonzero(rare_indices[k:] - np.arange(ri, ri + len(rare_indices[k:]), 1))[0]
                    if len(nnind) == 0:
                        nnind = np.nonzero(rare_indices[:k+1] - np.arange(ri - k, ri+1, 1))[0]
                        if len(nnind) == 0:
                            continue
                        else:
                            index = ri - nnind[-1]
                    else:
                        index = ri + nnind[0]
                    degrees[ri] = degrees[index]

        unique = np.unique(degrees, axis=0)
        degree_max_index = np.nonzero(max_values[:max_degree] > max_values[:max_degree].max() * 0.75)[0][-1]
        sort_data = [unique[:,i] for i in max_values[max_degree:].argsort() + max_degree] + [unique[:,i] for i in range(degree_max_index + 1)]
        ind = np.lexsort(sort_data)
        unique = unique[ind][::-1]
        img = np.equal(degrees[:,np.newaxis], unique).all(axis=2).T

        if filtering > 0:
            if logging:
                print('Filter result')
            filtered = []
            length = img.shape[1]
            for i in range(length):
                si_length = min(length, i+filtering) - max(0, i-filtering)
                si = img[:,max(0, i-filtering):min(length, i+filtering)].sum(axis=1)
                si[si < si_length * 0.5] = 0
                si[si > 0] = 1
                filtered.append(si)
            img = np.array(filtered).T
            zero_indices = (img.sum(axis=-1) > threshold)
            img = img[zero_indices]
            unique = np.array(unique)[zero_indices]

        if logging:
            print('Create labels')
        unique = unique.astype('<U8')
        unique = np.insert(unique, max_degree, np.array(['|'] * len(unique)), axis=1)
        zero_indices = ~np.all(unique == '0', axis=0)
        unique = unique[:, zero_indices]
        unique_labels = np.array(unique_labels)[zero_indices]
        unique = [' '.join(u) for u in unique]
        unique_labels = ' '.join(unique_labels)

        return img, unique, unique_labels


    """ Calculate the chain permutation invariant RMSD between the initial conformation and the structure at frame id `frame`
    The method tries different rotations of the structure to find the best chain permutation. However, this does not
    guarantee that the optimal permutation is found, especially for incomplete assemblies.

    Keyword arguments:
    frame -- index of frame in trajectory to use for the calculation
    angle_steps_per_axis -- number of steps per euler angle for the rotation trials (computational effort is n^3)
    overwrite_chain_slices -- set the amino acid selection slice (value) for chain at index (key)
    """
    def permutation_invariant_self_rmsd(self, frame: int, angle_steps_per_axis: int=30, overwrite_chain_slices: dict=None):
        slices = self._get_chain_slices(overwrite=overwrite_chain_slices)
        number_of_chains = self.get_number_of_chains()
        indentical_chains = np.ones((number_of_chains, number_of_chains), dtype=bool)
        for i in range(number_of_chains):
            for j in range(i+1, number_of_chains):
                if self._residue_names[slices[i]] != self._residue_names[slices[j]]:
                    indentical_chains[i,j] = indentical_chains[j,i] = False
        init_coordinates = np.empty((number_of_chains, 3))
        frame_coordinates = np.empty((number_of_chains, 3))
        for i in range(number_of_chains):
            init_coordinates[i] = self._initial_conformation[slices[i]].mean(axis=0)
            frame_coordinates[i] = self._chains[frame, slices[i]].mean(axis=0)
        init_center = init_coordinates.mean(axis=0)
        frame_center = frame_coordinates.mean(axis=0)
        init_coordinates_centered = init_coordinates - init_center
        frame_coordinates_centered = frame_coordinates - frame_center
        def coarse_min_rmsd(z_rot, y_rot, x_rot, degrees=True):
            rotation = Rotation.from_euler('zyx', [z_rot, y_rot, x_rot], degrees=degrees)
            rotated = rotation.apply(frame_coordinates_centered)
            distances = cdist(init_coordinates_centered, rotated, 'sqeuclidean')
            distances += (1 - indentical_chains) * distances.max()
            indices = distances.argmin(axis=1)
            """if z_rot + y_rot + x_rot == 0: # Debugging
                print(distances)
                print(indices)
                print(distances[len(indices), indices])"""
            distance_sum = distances[np.arange(len(indices)), indices].sum()
            return distance_sum, indices
        min_values = {}
        angle_steps = np.linspace(0, 360, angle_steps_per_axis, endpoint=False)
        for z,y,x in itertools.product(angle_steps, angle_steps, angle_steps):
            distance_sum, indices = coarse_min_rmsd(z, y, x, degrees=True)
            key = ' '.join(map(str,indices))
            if key in min_values:
                if distance_sum < min_values[key]:
                    min_values[key] = distance_sum
            else:
                min_values[key] = distance_sum
        min_info = sorted(list(min_values.items()), key=lambda e: e[1])
        min_rmsd = 10e10
        result_indices = []
        for min_i in min_info[:10]:
            min_indices = list(map(int, min_i[0].split()))
            permutated_frame_coordinates = []
            for i in min_indices:
                permutated_frame_coordinates.append(self._chains[frame, slices[i]])
            permutated_frame_coordinates = np.concatenate(permutated_frame_coordinates, axis=0)
            rmsd = aligned_rmsd(self._initial_conformation, permutated_frame_coordinates)
            if rmsd < min_rmsd:
                min_rmsd = rmsd
                result_indices = min_indices
        return min_rmsd, result_indices
        
def aligned_rmsd(coordinates1: np.ndarray, coordinates2: np.ndarray) -> float:
    """ Calculate the RMSD between two lists of coordinates after their translational and rotational alignment """
    assert coordinates1.shape == coordinates2.shape, 'Both input arrays must have the same shape'
    assert len(coordinates1.shape) == 2 and coordinates1.shape[-1] == 3, 'The input array size must be (N, 3)'
    center1 = coordinates1.mean(axis=0)
    center2 = coordinates2.mean(axis=0)
    centered1 = coordinates1 - center1
    centered2 = coordinates2 - center2
    rot, rssd = Rotation.align_vectors(centered1, centered2)
    rmsd = rssd / np.sqrt(len(coordinates1))
    return rmsd

def read_xvg_file(filename: str, dataframe: bool=True):
    """ Read xvg file as np.ndarray or pd.DataFrame """
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    labels = []
    for line in lines:
        if line[0] == '@' and 'xaxis' in line:
            labels.append(line.split('"')[1])
        elif line[0] == '@' and 'legend "' in line:
            labels.append(line.split('"')[1])
    data = np.array([[float(x) for x in line.split() if x] for line in lines if line[0] != '#' and line[0] != '@'])
    while len(labels) < data.shape[1]:
        labels.append(len(labels))
    if dataframe:
        return pd.DataFrame(data=data, columns=labels)
    else:
        return data
    
