from argparse import ArgumentParser
import glob
import multiprocessing as mp
import os
import lzma
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


aa_3c_to_1c = {
    "ALA":'A',
    "CYS":'C',
    "ASP":'D',
    "GLU":'E',
    "PHE":'F',
    "GLY":'G',
    "HIS":'H',
    "ILE":'I',
    "LYS":'K',
    "LEU":'L',
    "MET":'M',
    "ASN":'N',
    "PRO":'P',
    "GLN":'Q',
    "ARG":'R',
    "SER":'S',
    "THR":'T',
    "VAL":'V',
    "TRP":'W',
    "TYR":'Y',
}


def distribute(lst:list, n_bins:int) -> list:
    """
        Returns a list containg n_bins number of lists that contains the items passed in with the 1st argument

        :param lst: list that contains that items to be distributed across n bins
        :param n_bins: number of bins/lists across which to distribute the items of lst
    """ 
    if n_bins < 1:
       raise ValueError('The number of bins must be greater than 0')
    
    #cannot have empty bins so max number of bin is always less than or equal to list length
    n_bins = min(n_bins, len(lst))
    distributed_lists = []
    for i in range(0, n_bins):
        distributed_lists.append([])
    
    for i, item in enumerate(lst):
        distributed_lists[i%n_bins].append(item)

    return distributed_lists


def atom_from_pdb_line(atom_line:str) -> dict:
    """
        Parses a single line string in the standard PDB format and returns a list a dict that represents an atom with 3d coordinates and a type(element)

        :param atom_line: string representing a single line entry in a PDB file which contains information about an atom in a protein structure
    """ 
    coordinates = np.array([float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54])])
    return {"type":atom_line[13:16].strip(),"xyz":coordinates,}


def get_lines_from_pdb_file(pdb_filename:str) -> list:
    """
        Returns the contents of a protein databank file (PDB) file as a list of strings

        :param pdb_filename: string representing the path of the PDB file to open and parse (can handle PDB files that have been compressed via GZIP or LZMA)
    """ 

    if not os.path.isfile(pdb_filename):
        raise ValueError('Non existing PDB file was specified')

    pdb_file = None 
    if(pdb_filename.endswith('.xz')):
        pdb_file = lzma.open(pdb_filename, 'rt')
    elif(pdb_filename.endswith('.gz')):
        pdb_file = gzip.open(pdb_filename,'rt')
    elif(pdb_filename.endswith('.pdb')):
        pdb_file = open(pdb_filename,'rt')
    else:
        raise ValueError('Unable to parse a PDB file with invalid file extension')
    
    pdb_data = pdb_file.read()
    pdb_file.close()
    return pdb_data.splitlines()



def get_distogram_data(pdb_filepath):

    coords = []
    chain_ranges = []
    active_range = None
    residues = []
    prev_chain = None

    res_count = 0

    for atom_line in get_lines_from_pdb_file(pdb_filepath):
        if atom_line[0:4] != 'ATOM' or atom_line[13:16].strip() != 'CA' :
            continue

        res_count += 1

        chain = atom_line[20:22].strip()
        if prev_chain != chain:
            active_range = {'start':res_count, 'end':res_count}
            chain_ranges.append(active_range)

        prev_chain = chain

        active_range['end'] = res_count

        aa_type = aa_3c_to_1c[atom_line[17:20]]
        residues.append(aa_type)

        atom = atom_from_pdb_line(atom_line)
        coords.append(atom['xyz'])

    num_ca_atoms = len(coords)
    c1_matrix = np.tile(coords, (1, num_ca_atoms)).reshape(num_ca_atoms, num_ca_atoms, 3)
    c2_matrix = np.tile(coords, (num_ca_atoms, 1)).reshape(num_ca_atoms, num_ca_atoms, 3)
    distances = np.sqrt(np.sum((c1_matrix - c2_matrix)**2, axis=2))
    return distances, chain_ranges, residues

def save_heatmap_image_for_distogram(d, save_path, lines=[]):
   
    plt.figure(figsize=(8, 6))
    heatmap = plt.imshow(d, cmap='YlGnBu_r', vmin=0, vmax=30, interpolation='nearest')
    # Add color bar at the side
    plt.colorbar(label='CA Angstroms')
    for pos in ['right', 'top', 'bottom', 'left']: 
        plt.gca().spines[pos].set_visible(False) 
    # Add labels for clarity
    plt.xlabel('Residue 1')
    plt.ylabel('Residue 2')

    for x in lines:
        plt.axvline(x=x, color='black', linestyle='-')

    for y in lines:
        plt.axhline(y=y, color='black', linestyle='-')

    # Add title
    plt.title('CA distogram')
    plt.savefig(save_path)
    plt.close()



def compute_distograms(output_folder, cpu_index, pdb_files):

    for pdb_filepath in pdb_files:

        print(f"starting {pdb_filepath}")
        distogram, chain_ranges, residues  = get_distogram_data(pdb_filepath)
        
        lines = []
        for i in range(0, len(chain_ranges) -1):
            r1 = chain_ranges[i]
            r2 = chain_ranges[i + 1]
            lines.append(0.5*(r1['end'] + r2['start']))

        filename = pdb_filepath.split('/').pop()

        distogram_df = pd.DataFrame(distogram).round(1)
        new_row_df = pd.DataFrame([residues], columns=distogram_df.columns)
        distogram_df = pd.concat([new_row_df, distogram_df]).reset_index(drop=True)

        distogram_df.insert(0, 'residues', [' '] + [residue for residue in residues])
        distogram_df.to_csv(os.path.join(output_folder, f"{filename}.pdb_distogram.csv"), index=None)
        save_heatmap_image_for_distogram(distogram, os.path.join(output_folder, f"{filename}.pdb_distogram.pdf"), lines)
        print(f"finished {pdb_filepath}")


def analysis_thread_did_finish(arg1):
   return

def analyze_folder(data_folder:str, name_filter:str,) -> str:

    data_folder = data_folder.rstrip('/')
    output_folder = os.path.basename(data_folder) + "_distograms"
    
    wildchar_str = '*.pdb.*'
    if name_filter is not None and len(name_filter) > 0:
        wildchar_str = '*' + name_filter + wildchar_str

    pdb_filepaths = glob.glob(os.path.join(data_folder,wildchar_str))
    if len(pdb_filepaths) < 1:
        print(f"FOLDER {data_folder} does not contain any PDB files matching the filer criteria")
        return None


    index = 1
    while os.path.isdir(output_folder):
        #if we find existing folders with the output folder name we will iterate over index until we find an unused folder name
        output_folder = os.path.basename(data_folder) + "_distograms_" + str(index)
        index += 1

    #guaranteed to have a new unique output_folder name, lets make it     
    os.mkdir(output_folder)

    num_cpus_to_use = mp.cpu_count()
    num_cpus_to_use = min(num_cpus_to_use,len(pdb_filepaths))
    print(f"Splitting analysis job across {num_cpus_to_use} different CPUs")
    pool = mp.Pool(num_cpus_to_use)
    pdb_filepath_lists = distribute(pdb_filepaths, num_cpus_to_use)

    traces = []
    for cpu_index in range(0, num_cpus_to_use):
        #create a new thread to analyze the complexes (1 thread per CPU)
        traces.append(pool.apply_async(compute_distograms, 
                                       args=(output_folder, cpu_index, pdb_filepath_lists[cpu_index]), 
                                       callback=analysis_thread_did_finish)
        )

    for t in traces:
        t.get()

    pool.close()
    pool.join()

    return output_folder



if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        "input",
        default="",
        help="One or more folders with PDB files(can be compressed in GZIP OR LZMA format).",
        nargs='*',)    
    parser.add_argument(
        "--name-filter",
        default='',
        help="Filter structures by name so that only files containing this string/name are processed",
        type=str,)
    

    output_folders = []
    args = parser.parse_args()
    for folder in args.input:
    
        if not os.path.isdir(folder):
            print(f"ERROR {folder} does not appear to be a non valid folder, skipping")
            continue

        print(f"Starting to analyze folder ({folder})")
        output_folders.append(analyze_folder(folder, args.name_filter))
        print(f"Finished analyzing folder ({folder})")
        print(" "*80)
        print("*"*80)
        print("*"*80)
        print(" "*80)