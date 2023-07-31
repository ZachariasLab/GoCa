from re import S
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from scipy.special import expit
from scipy import optimize
import networkx as nx
from process_helper import parmap
import matplotlib.pyplot as plt

class Trajectory:

    def __init__(self):
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

    def __len__(self):
        return len(self._chains)

    def size(self, chainId=None):
        start_index, stop_index = self._get_chain_start_stop_indices(chainId)
        return stop_index - start_index

    def get_number_of_chains(self):
        return len(self._chain_indices)

    def get_time_steps(self):
        return self._time_steps

    def get_residue_names(self):
        return self._residue_names

    def get_box_dimensions(self):
        return self._box_dimensions

    def get_contacts_inside_chains(self):
        return self._contacts_inside_chains

    def get_contacts_between_chains(self):
        return self._contacts_between_chains

    def _get_chain_start_stop_indices(self, chainId):
        if chainId == None:
            return 0, self._chains.shape[1]
        start_index = self._chain_indices[chainId]
        stop_index = self._chains.shape[1] if (chainId == len(self._chain_indices) - 1) else self._chain_indices[chainId + 1]
        return start_index, stop_index

    def __getitem__(self, pos):
        if type(pos) is int:
            start_index, stop_index = self._get_chain_start_stop_indices(pos)
            return self._chains[:, start_index:stop_index]
        elif type(pos) is tuple and len(pos) == 2:
            start_index, stop_index = self._get_chain_start_stop_indices(pos[0])
            return self._chains[pos[1], start_index:stop_index]
        elif type(pos) is tuple and len(pos) == 3:
            start_index, stop_index = self._get_chain_start_stop_indices(pos[0])
            return self._chains[pos[1], start_index:stop_index][pos[2]]

    def is_valid(self):
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

    def _get_single_fraction_of_native_contact(self, coordinates, native_contacts, cut_off, sigmoid_factor, pbc):
        a1 = coordinates[native_contacts[:,0]]
        a2 = coordinates[native_contacts[:,1]]
        if pbc:
            delta = np.abs(a1 - a2)
            distances = np.linalg.norm(np.where(delta > 0.5 * self._box_dimensions, delta - self._box_dimensions, delta), axis=1)
        else:
            distances = np.linalg.norm(a1 - a2, axis=1)
        if sigmoid_factor <= 0:
            return np.mean(distances < cut_off)
        else:
            return np.mean(expit(sigmoid_factor*(cut_off-distances)))
    
    def _get_fraction_of_native_contacts(self, coordinates, native_contacts, cut_off, sigmoid_factor, pbc):
        if coordinates.shape[0] == 1:
            return self._get_single_fraction_of_native_contact(coordinates[0], native_contacts, cut_off, sigmoid_factor, pbc)
        else:
            result = []
            result = np.array(parmap(
                lambda frame: self._get_single_fraction_of_native_contact(frame, native_contacts, cut_off, sigmoid_factor, pbc), coordinates))
            return np.array(result)

    def get_fraction_of_native_contacts_inside(self, frame=None, chain_id=None, offset=0.15, sigmoid_factor=20, pbc=False):
        frame_selection = Trajectory.__frame_to_slice(frame)
        chain_start, chain_stop = self._get_chain_start_stop_indices(chain_id)
        selector = (self._contacts_inside_chains >= chain_start).all(axis=1) & (self._contacts_inside_chains < chain_stop).all(axis=1)
        contacts = self._contacts_inside_chains[selector]
        cut_off = self._contacts_inside_chains_distances[selector] + offset
        return self._get_fraction_of_native_contacts(self._chains[frame_selection], contacts, cut_off, sigmoid_factor, pbc)

    def get_fraction_of_native_contacts_between(self, frame=None, offset=0.15, sigmoid_factor=20, pbc=False):
        frame_selection = Trajectory.__frame_to_slice(frame)
        digitized = np.digitize(self._contacts_between_chains, self._chain_indices) - 1
        contacts = self._contacts_between_chains[digitized[:,0] != digitized[:,1]]
        cut_off = self._contacts_between_chains_distances[digitized[:,0] != digitized[:,1]] + offset
        return self._get_fraction_of_native_contacts(self._chains[frame_selection], contacts, cut_off, sigmoid_factor, pbc)

    def get_fraction_of_native_contacts_chain_pairs(self, frame=None, offset=0.15, sigmoid_factor=20, pbc=False):
        frame_selection = Trajectory.__frame_to_slice(frame)
        digitized = np.digitize(self._contacts_between_chains, self._chain_indices) - 1
        groups = np.unique(digitized, axis=0)
        result = []
        for group in groups:
            contacts = self._contacts_between_chains[(digitized[:,0] == group[0]) & (digitized[:,1] == group[1])]
            cut_off = self._contacts_inside_chains_distances[(digitized[:,0] == group[0]) & (digitized[:,1] == group[1])] + offset
            result.append(self._get_fraction_of_native_contacts(self._chains[frame_selection], contacts, cut_off, sigmoid_factor, pbc))
        return result, groups

    def _get_chain_slices(self, overwrite=None):
        indices = [*self._chain_indices, len(self._residue_names)]
        slices = [slice(indices[i], indices[i+1], 1) for i in range(len(self._chain_indices))]
        if overwrite is not None and type(overwrite) == dict:
            for key, value in overwrite.items():
                slices[key] = slice(slices[key].start + value.start, slices[key].start + value.stop, value.step if value.step else 1)
        return slices

    def get_chain_distance_map(self, pbc=False, overwrite_chain_slices=None):
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

    def get_chain_graph_helper(self, rolling_window_size=40, overwrite_chain_slices=None, processes=None):
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

    def get_chain_graph(self, rolling_window_size=30, std_cut_off=0.3, distance_cutoff_factor=0.5, single_distance_cutoff=1.2,
            fixed_distance_cutoff=None, overwrite_chain_slices=None, logging=True, show_distance_histogram=False, processes=None):
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

        
    def get_chain_formations(self, graphs, threshold=100, filtering=0, logging=False):
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
        

def aligned_rmsd(coordinates1, coordinates2):
    center1 = coordinates1.mean(axis=0)
    center2 = coordinates2.mean(axis=0)
    centered1 = coordinates1 - center1
    centered2 = coordinates2 - center2
    rot, rmsd = Rotation.align_vectors(centered1, centered2)
    return rmsd

def get_fitted_unfolding_curve(temperatures, rmsd, p0=None, error=False, bounds=None, maxfev=1000, method=None):
    if p0 is None:
        p0 = [max(rmsd), 0.2, np.median(temperatures), min(rmsd)]
    np.seterr(over='ignore')
    def sigmoid(x, a, b, c, d):
        return a / (1 + np.exp(-b*(x-c))) + d
    if bounds is None:
        params, cov = optimize.curve_fit(sigmoid, temperatures, rmsd, p0=p0, maxfev=maxfev, method=method)
    else:
        params, cov = optimize.curve_fit(sigmoid, temperatures, rmsd, p0=p0, bounds=bounds, maxfev=maxfev, method=method)
    tpx = params[2]
    np.seterr(over='warn')
    if error:
        return tpx, lambda x: sigmoid(x, *params), np.sqrt(cov[2,2])
    else:
        return tpx, lambda x: sigmoid(x, *params)

def read_xvg_file(filename, dataframe=True):
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
    
