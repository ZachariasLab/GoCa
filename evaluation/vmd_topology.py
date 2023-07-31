"""
Run the following script (with the correct filepaths) or add it to your .vmdrc file:

proc loadSBCGMcontacts {file} {
  set result [exec python "/path_to_this_script/vmd_topology.py" $file]
  eval $result
}
loadSBCGMcontacts "/path_to_the_topology_file_of_the_analysed_structure.top"
"""

import sys
import os.path
import numpy as np

def main():
    if len(sys.argv) < 2 or len(sys.argv[1]) < 1:
        print('puts "Error: Please specify a topology file"')
        exit()

    filename = sys.argv[1]
    print('puts "Loading: ' + filename + '"')

    if not os.path.isfile(filename):
        print('puts "Error: File does not exist"')
        exit()

    file = open(filename, 'r')

    linestack = []
    status = 0
    non_bonded = []
    molecules = []
    atom_indices = {}
    atom_types = {}
    name = ''
    inter_pairs = []
    atom_types_all = []
    chain_ids = []
    chain_count = 0
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
            elif '[ molecules ]' in data:
                status = 3
            elif '[ moleculetype ]' in data:
                status = 4
            elif len(data) < 2:
                status = 0
            elif status == 1:
                non_bonded.append([x for x in data.split() if x][:2])
            elif status == 4:
                name = data.split()[0]
                molecules.append(name)
                atom_indices[name] = []
                atom_types[name] = []
            elif status == 2:
                data_split = [x for x in data.split() if x]
                atom_indices[name].append(int(data_split[0]))
                atom_types[name].append(data_split[1])
            elif status == 3:
                data_split = [x for x in data.split() if x]
                m_name = data_split[0]
                number = int(data_split[1])
                for m in range(number):
                    atom_types_all.extend(atom_types[m_name])
                    chain_ids.extend(list(np.ones(len(atom_types[m_name])) * chain_count))
                    chain_count += 1
    atom_types_all = np.array(atom_types_all)
    for nb in non_bonded:
        for i1 in np.where(atom_types_all == nb[0])[0]:
            for i2 in np.where(atom_types_all == nb[1])[0]:
                if chain_ids[i1] != chain_ids[i2]:
                    inter_pairs.append([i1, i2])
    file.close()

    print('label delete Bonds all')
    for pair in inter_pairs:
        print('set distance [measure bond {{{} {}}}]'.format(*pair))
        print('if {$distance < 10} {\nlabel add Bonds ' + '0/{} 0/{}'.format(*pair) + '\n}')


if __name__ == "__main__":
    main()