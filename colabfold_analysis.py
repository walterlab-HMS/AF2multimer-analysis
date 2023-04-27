#!/usr/bin/python3

from argparse import ArgumentParser
import glob
import multiprocessing as mp
import os
import re
import lzma
import gzip
import numpy as np
import math
from statistics import mean
import pandas as pd
import sys
import shutil

#dict for converting 3 letter amino acid code to 1 letter code
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

def join_csv_files(files:list, output_name:str, sort_col:str = None, sort_ascending:bool = False, headers = None):
    """
        Join multiple CSV files into a single file.

        :param files (list): A list of file paths to CSV files to be joined.
        :param output_name (str): The name of the output file.
        :param headers (list, optional): A list of column names for the output file. If not provided, the column names from the first input file are used.
    """
    if(len(files) < 1):
        return

    first_df = pd.read_csv(files[0])
    combo_df = pd.DataFrame(columns = first_df.columns)
    for f in files:
        df_temp = pd.read_csv(f)
        combo_df = pd.concat([combo_df, df_temp], ignore_index=True)

    if headers is not None:
        combo_df.columns = headers

    if sort_col:
        combo_df.sort_values(by=[sort_col], ascending=sort_ascending, inplace=True)
    combo_df.to_csv(output_name, index=None)


def distribute(lst:list, n_bins:int) -> list:
    """
        Returns a list containg n_bins number of lists that contains the items passed in with the lst argument

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


def get_af_model_num(filename) -> int:
    """
        Returns the Alphafold model number from an input filestring as an int

        :param filename: string representing the filename from which to extract the model number
    """ 
    
    if "model_" not in filename: return 0

    model_num = int(re.findall(r'model_\d+', filename)[0].replace("model_", ''))
    return model_num


def get_finished_complexes(path:str, search_str:str = '.done.txt') -> list:
    """
        Returns a list of string values representing the names of the complexes that were found in a specified path (folder)

        :param path: string representing the path/folder to search for complexes
        :param search_str: string representing the pattern to use when searching for the completed complexes. By default this is *.done.txt because Colabfold outputs 1 such file per complex.
    """

    done_files = glob.glob(os.path.join(path,'*' + search_str))
    complex_names = [os.path.basename(f).replace(search_str, '') for f in done_files]
    return complex_names


def get_filepaths_for_complex(path:str, complex_name:str, pattern:str = '*') -> list:
    """
        Helper methdof for returning a list of filepaths (strs) that match the specified GLOB pattern

        :param path: string representing the path/folder to search for complexes
        :param complex_name: string that represents the name of a complex
        :param pattern: string representing the pattern to use when searching for files belonging to complex. Ex: *.pdb, *.json, etc
    """

    glob_str = os.path.join(path, complex_name + pattern)
    return sorted(glob.glob(glob_str))


def get_pae_values_from_json_file(json_filename) -> list:
    """
        Returns a list of string values representing the pAE(predicated Aligned Error) values stored in the JSON output

        :param json_filename: string representing the JSON filename from which to extract the PAE values
    """ 

    if not os.path.isfile(json_filename):
        raise ValueError('Non existing PAE file was specified')

    scores_file = None 
    if(json_filename.endswith('.xz')):
        scores_file = lzma.open(json_filename, 'rt')
    elif(json_filename.endswith('.gz')):
        scores_file = gzip.open(json_filename,'rt')
    elif (json_filename.endswith('.json')):
        scores_file = open(json_filename,'rt')
    else:
        raise ValueError('pAE file with invalid extension cannot be analyzed. Only valid JSON files can be analyzed.')

    #read pae file in as text
    file_text = scores_file.read()
    pae_index = file_text.find('"pae":')
    scores_file.close()
    
    #Transform string representing 2d array into a 1d array of strings (each string is 1 pAE value). We save time by not unecessarily converting them to numbers before we use them.
    pae_data = file_text[pae_index + 6:file_text.find(']]', pae_index) + 2].replace('[','').replace(']','').split(',')

    if len(pae_data) != int(math.sqrt(len(pae_data)))**2:
        #all valid pAE files consist of an N x N matrice of scores
        raise ValueError('pAE values could not be parsed from files')
    
    return pae_data

def get_pdockq_val(num_contacts:int, avg_plddt:float) -> float:
    """
        Returns the predicted DOCKQ value, a parameter ranging from 0 to 1 that estimates how well a predicted interfaces matches a true interface

        :param num_contacts: the number of contacts(unique residue pairs) in the interface
        :param avg_plddt: the average pLDDT value over the interface
    """ 
    #formula represents a sigmoid function that was empirically derived in the paper here: https://www.nature.com/articles/s41467-022-28865-w
    #even though pDOCKQ range is supposed to bve from 0 to 1, this function maxes out at 0.742
    return 0.724/(1 + 2.718**(-0.052*(avg_plddt*math.log(num_contacts, 10) - 152.611))) + 0.018


def dist2(v1, v2) -> float:
    """
        Returns the square of the Euclian distance between 2 vectors carrying 3 values representing positions in the X,Y,Z axis

        :param v1: a vector containing 3 numeric values represening X, Y, and Z coordinates
        :param v2: a vector containing 3 numeric values represening X, Y, and Z coordinates
    """ 

    if len(v1) != 3:
        raise ValueError('3D coordinates require 3 values')
    
    if len(v2) != 3:
        raise ValueError('3D coordinates require 3 values')

    return (v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2

def atom_from_pdb_line(atom_line:str) -> dict:
    """
        Parses a single line string in the standard PDB format and returns a list a dict that represents an atom with 3d coordinates and a type(element)

        :param atom_line: string representing a single line entry in a PDB file which contains information about an atom in a protein structure
    """ 
    coordinates = np.array([float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54])])
    return {"type":atom_line[13:16].strip(),"xyz":coordinates,}

def get_closest_atoms(res1:dict, res2:dict):
    """
        Find the two closest atoms between two residues and returns the minimum distance as well as the closest atoms from each residue.

        :param res1: A dictionary representing the first residue. It should have a key 'atoms' that contains a list of dictionaries, where each dictionary represents an atom with keys 'xyz' (a list of three floats representing the x, y, and z coordinates) and 'name' (a string representing the name of the atom).
        :param res2: A dictionary representing the second residue. It should have a key 'atoms' that contains a list of dictionaries, where each dictionary represents an atom with keys 'xyz' (a list of three floats representing the x, y, and z coordinates) and 'name' (a string representing the name of the atom).
    """
    min_d2 = 1e6
    atoms = [None, None]
    for a1 in res1['atoms']:
        for a2 in res2['atoms']:
            d2 = dist2(a1["xyz"], a2["xyz"])
            if d2 < min_d2:
                min_d2 = d2
                atoms[0] = a1
                atoms[1] = a2
    return (min_d2, atoms)

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

def get_contacts_from_structure(pdb_filename:str, max_distance:float = 8, min_plddt:float = 70) -> dict:
    """
        Returns a dict that contains all amino acids between different chains that are in contact and meet the specified criteria

        :param pdb_filename:string of the filepath to the PDB structure file to be parsed for contacts
        :param max_distance:the maximum allowed distance in Angstroms that the 2 residues must have in order to be considered in contact
        :param min_plddt:the minimum pLDDT(0-100) 2 residues must have to be considered in contact
    """ 

    #holds a more relaxed distance criteria for fast preliminary filtering  
    d2_n_cutoff = (max_distance + 20)**2

    d2_cutoff = max_distance**2
    last_chain = None
    abs_res_index = 0
    chain_index = -1
    chains = []

    residues = []

    #holds 3d coordinates of all amide nitrogens for all residues
    #organized as a 2d list with rows representing chains and columns are residues within the chain
    N_coords = []
    
    for atom_line in get_lines_from_pdb_file(pdb_filename):
        if atom_line[0:4] != 'ATOM':
            continue
        
        atom_type = atom_line[13:16].strip()
        is_nitrogen = atom_type == 'N'
        if(is_nitrogen):
            #Keep track of what the absolute index of this residue in the file is 
            abs_res_index += 1

        #in AlphaFold output PDB files, pLDDT values are stored in the "bfactor" column 
        bfactor = float(atom_line[60:66])
        if bfactor < min_plddt:
            #No need to examine this atom since it does not meet the pLDDT cutoff, skip it
            continue

        atom = atom_from_pdb_line(atom_line)
        if is_nitrogen:
            #Every amino acid residue starts PDB entry with exactly one "N" atom, so when we see one we know we have just encountered a new residue
            chain = atom_line[20:22].strip()
            if chain != last_chain:
                #Are we in a new chain? If so, increment chain index and create new list in "residuesa"
                chain_index += 1
                last_chain = chain
                N_coords.append([])
                residues.append([])
                chains.append(chain)

            residue = {"chain":chain, "atoms":[],'c_ix':int(atom_line[22:26]), "a_ix":abs_res_index, "type":aa_3c_to_1c[atom_line[17:20]], "plddt":bfactor}
            residues[chain_index].append(residue)

            #add nitrogen atom coordinates to coordinates list to allow for fast broad searching later
            N_coords[chain_index].append(atom['xyz'])

        residue['atoms'].append(atom)
    
    contacts = []
    num_chains = len(chains)

    #loop through all the protein chains to find contacts between chains
    #strategy is to first look for residues in general proximity by just looking at the distance between their amide nitrogen
    for i in range(0, num_chains):
        chain_1_coords = N_coords[i]
        num_in_c1 = len(chain_1_coords)
        for i2 in range(i + 1, num_chains):

            chain_2_coords = N_coords[i2]
            num_in_c2 = len(chain_2_coords)

            #construct 2 3D numpy arrays to hold coordinates of all residue amide nitorgens
            c1_matrix = np.tile(chain_1_coords, (1, num_in_c2)).reshape(num_in_c1, num_in_c2, 3)
            c2_matrix = np.tile(chain_2_coords, (num_in_c1, 1)).reshape(num_in_c1, num_in_c2, 3)

            #calculate euclidian distance squared (faster) between all amide nitorgens of all residues
            d2s = np.sum((c1_matrix - c2_matrix)**2, axis=2)
            #get residue pairs where amide nitrogens are closer than the initial broad cutoff
            index_pairs = list(zip(*np.where(d2s < d2_n_cutoff)))

            #find closest atoms between residues that were found to be somewhat in proximity
            for c1_res_ix, c2_res_ix in index_pairs:
                
                r1 = residues[i][c1_res_ix]
                r2 = residues[i2][c2_res_ix]
                min_d2, atoms = get_closest_atoms(r1, r2)
                if(min_d2 < d2_cutoff):
                    #residues have atoms closer than specified cutoff, lets add them to the list
                    contacts.append({
                        'distance': round(math.sqrt(min_d2), 1),
                        "aa1":{"chain":r1["chain"], "type":r1["type"], "c_ix":r1['c_ix'], "a_ix":r1['a_ix'], "atom":atoms[0]['type'], "plddt": r1["plddt"]},
                        "aa2":{"chain":r2["chain"], "type":r2["type"], "c_ix":r2['c_ix'], "a_ix":r2['a_ix'], "atom":atoms[1]['type'], "plddt": r2["plddt"]}
                    })
    return contacts


def get_contacts(pdb_filename:str, pae_filename:str, max_distance:float, min_plddt:float, max_pae:float, pae_mode:str) -> dict:
    """
        Get contacts from a protein structure in PDB format that meet the specified distance and confidence criteria.

        :param pdb_filename (str): The path to the PDB file.
        :param pae_filename (str): The path to the predicted Alignment Error (pAE) file.
        :param max_distance (float): The maximum distance between two atoms for them to be considered in contact.
        :param min_plddt (float): The minimum PLDDT score required for a residue to be considered "well-modeled".
        :param max_pae (float): The maximum predicted Alignment Error allowed for a residue to be considered "well-modeled".
        :param pae_mode (str): The method to use for calculating predicted atomic error (pAE). Possible values are "avg" or "min".
    """

    model_num = get_af_model_num(pdb_filename)
    if model_num < 1 or model_num > 5:
        raise ValueError('There are only 5 Alphafold models, numbered 1 to 5. All PDB files and PAE files must have a valid model number to be analyzed.')

    if model_num != get_af_model_num(pae_filename):
        raise ValueError('File mismatch, can only compare PDB and PAE files from same complex and the same AF2 model')

    #first determine which residues are in physical contact(distance) and have a minimum pLDDT score (bfactor column)
    contacts = get_contacts_from_structure(pdb_filename, max_distance, min_plddt)
    if len(contacts) < 1:
        return {}
    
    #extract PAE data as a list of strings("PAE values") from the PAE file which is a linearized form of a N by N matrix where N is the total number of residues inb the predicted structure 
    pae_data = get_pae_values_from_json_file(pae_filename)

    #need this value for converting between amino acid index and the pAE array index
    total_aa_length = int(math.sqrt(len(pae_data)))

    filtered_contacts = {}

    for c in contacts:

        aas = [c['aa1'], c['aa2']]
        aa_indices = [aas[0]['a_ix'],aas[1]['a_ix']]

        #convert the absolute amino acid index into the linear index where the 2 PAE values for each pair are (x, y) and (y, x)
        pae_index_1 = total_aa_length*(aa_indices[0] - 1) + aa_indices[1] - 1
        pae_index_2 = total_aa_length*(aa_indices[1] - 1) + aa_indices[0] - 1

        #pae data contains string values, have to convert them to floats before using them for math calculations
        pae_values = [float(pae_data[pae_index_1]), float(pae_data[pae_index_2])]

        pae_value = 0.5*(pae_values[0] + pae_values[1]) if pae_mode == 'avg' else min(pae_values[0],pae_values[1])
        
        if(pae_value > max_pae):
            #The pAE value of this residue pair is too high, skip it
            continue

        #This contact meets all the specified criteria, add it to the dict

        #Use the 2 chains IDS as a key
        chain_contact_id = aas[0]['chain'] + ":" + aas[1]['chain']
        if chain_contact_id not in filtered_contacts:
            filtered_contacts[chain_contact_id] = {}

        #Use the absolute indices of the two residues in the PDB file as the unique key for this pair/contact
        contact_id = str(aa_indices[0]) + '&' + str(aa_indices[1])
        filtered_contacts[chain_contact_id][contact_id] = {

            'chains':[aas[0]['chain'], aas[1]['chain']],
            'inchain_indices':[aas[0]['c_ix'], aas[1]['c_ix']],
            'types':[aas[0]['type'], aas[1]['type']], 
            'pae':pae_value,
            'paes':pae_values,
            'plddts':[aas[0]['plddt'], aas[1]['plddt']], 
            'model':model_num,
            'distance':c['distance']
        }

    return filtered_contacts

        
def calculate_interface_statistics(contacts:dict) -> dict:
    """
        Returns summary confidence statistics such as pAE and pLDDT values across all the contacts in an interface

        :param contacts:dict of contacts in an interface of the form {'chain1:chain2':{'1&400':{plddts:[75,70], paes:[10, 7]}, '4&600':{plddts:[68,77], paes:[8, 3]}}}
    """ 

    #plddts always range from 0 to 100
    plddt_sum = 0 
    plddt_min = 100
    plddt_max = 0
    plddt_avg = 0
   
    #paes always range from 0 to 30
    pae_avg = 0
    pae_sum = 0
    pae_min = 30
    pae_max = 0
    distance_avg = 0

    num_contacts = 0
    d_sum = 0

    for interchain_id, interchain_contacts in contacts.items():
        for contact_id, contact in interchain_contacts.items():

            avg_plddt = mean(contact['plddts'])
            plddt_sum += avg_plddt
            plddt_max = max(plddt_max, avg_plddt)
            plddt_min = min(plddt_min, avg_plddt)

            d_sum += contact['distance']

            pae_max = max(pae_max, contact['pae'])
            pae_min = min(pae_min, contact['pae'])
            pae_sum += contact['pae']

            num_contacts += 1

    pdockq = 0

    if num_contacts > 0:
        plddt_avg = round(plddt_sum/num_contacts, 1)
        pae_avg = round(pae_sum/num_contacts, 1)
        pdockq = round(get_pdockq_val(plddt_avg, num_contacts), 3)
        distance_avg = round(d_sum/num_contacts, 1)
    else:
        pae_min = 0
        plddt_min = 0

    data = {'pdockq':pdockq, 
            'num_contacts':num_contacts,
            'plddt':[plddt_min, plddt_avg, plddt_max],
            'pae':[pae_min, pae_avg, pae_max],
            'distance_avg': distance_avg}

    return data

def summarize_interface_statistics(interfaces:dict) -> dict:
    """
        summarize_interface_statistics returns aggregate statistics over multiple interfaces across predictions from different models

        :param interfaces:dict of interfaces in the form 
            {1: {'chain1:chain2':{
                                '1&400':{'plddts':[75,70], 'paes':[10, 7]}, 
                                '4&600':{'plddts':[68,77], 'paes':[8, 3]}
                                }
                }, 
            4: {'chain1:chain2':{
                                '13&400':{'plddts':[77,91], 'paes':[5, 7]}, 
                                '49&600':{'plddts':[68,56], 'paes':[9, 3]}
                                }
                },     
            }
    """
    
    unique_contacts = {}
    max_num_models = 0

    for model_num, interchain_interfaces in interfaces.items():
        for interchain_str, contacts in interchain_interfaces.items():
            for contact_id, c in contacts.items():

                if contact_id not in unique_contacts:
                    unique_contacts[contact_id] = 1
                else:
                    unique_contacts[contact_id] += 1

                max_num_models = max(max_num_models, unique_contacts[contact_id])

    num_contacts = 0
    sum_num_models = 0
    num_contacts_with_max_n_models = 0
    for contact_id, observation_count in unique_contacts.items():

        num_contacts += 1
        sum_num_models += observation_count

        if observation_count == max_num_models:
            num_contacts_with_max_n_models += 1
    
    summary_stats = {
        'max_n_models':max_num_models,
        'avg_n_models':round(sum_num_models/num_contacts, 1) if num_contacts > 0 else 0,
        'num_contacts_with_max_n_models':num_contacts_with_max_n_models,
        'num_unique_contacts':num_contacts
    }
    return summary_stats


def analyze_complexes(cpu_index:int, input_folder:str, output_folder:str, complexes:list, max_distance:float, min_plddt:float, max_pae:float, pae_mode:str):
    """
        Analyze protein complexes in PDB format.

        :param cpu_index (int): The index of the CPU being  used for parallel processing.
        :param input_folder (str): The path to the input folder containing PDB files.
        :param output_folder (str): The path to the output folder where analysis results will be saved.
        :param complexes (list): A list of complex names to be analyzed.
        :param max_distance (float): The maximum distance between two atoms for them to be considered in contact.
        :param min_plddt (float): The minimum PLDDT score required for a residue to be considered "well-modeled".
        :param max_pae (float): The maximum predicted alignment error allowed for a residue to be considered "well-modeled".
        :param pae_mode (str): The method to use for calculating predicted alignment error (PAE). Possible values are "min" or "avg".
    """

    summary_stats = {}
    interface_contacts = {}
    all_interface_stats = []
    all_contacts = []
    
    for index, cname in enumerate(complexes):

        print(f"Analyzing {index + 1} / {len(complexes)}: {cname}")

        pdb_filepaths = get_filepaths_for_complex(input_folder, cname, '*.pdb*')
        
        #get pAE files ending in .json or .json followed by two letters as would be the case of compressed gzipped files
        pae_filepaths = get_filepaths_for_complex(input_folder, cname, '*score*.json') + get_filepaths_for_complex(input_folder, cname, "*score*.json.??")
       
        if len(pdb_filepaths) != len(pae_filepaths):
            print(f"ERROR: Number of PDB files ({len(pdb_filepaths)}) does not match number of PAE files ({len(pae_filepaths)})")
            print("SKIPPING: " + cname)
            continue

        #sort the files by model number so that they are aligned for analysis ie PDB model 1 = PAE model 1
        pdb_filepaths.sort(key=get_af_model_num)
        pae_filepaths.sort(key=get_af_model_num)

        interface_contacts = {}

        #record which interface has the best score (ie the most high confidence contacts so we can report it out later)
        best_interface_stats = None

        for pdb_filename, pae_filename in zip(pdb_filepaths, pae_filepaths):
            
            model_num = get_af_model_num(pdb_filename)
            contacts = get_contacts(pdb_filename, pae_filename, max_distance, min_plddt, max_pae, pae_mode)
            interface_contacts[model_num] = contacts

            for interchain_str, interchain_interfaces in contacts.items():
                for contact_id, c in interchain_interfaces.items():
                    all_contacts.append({
                        "complex_name":cname,
                        "model_num":model_num,
                        "aa1_chain":c['chains'][0],
                        "aa1_index":c['inchain_indices'][0],
                        "aa1_type":c['types'][0],
                        "aa1_plddt":round(c['plddts'][0]),
                        "aa2_chain":c['chains'][1],
                        "aa2_index":c['inchain_indices'][1],
                        "aa2_type":c['types'][1],
                        "aa2_plddt":round(c['plddts'][1]),
                        "pae":c['pae'],
                    })

            if_stats = calculate_interface_statistics(contacts)

            if best_interface_stats is None:
                best_interface_stats = if_stats
                best_interface_stats['model_num'] = model_num
            else:

                if if_stats['pdockq'] > best_interface_stats['pdockq']:
                    best_interface_stats = if_stats
                    best_interface_stats['model_num'] = model_num

            all_interface_stats.append({
                "complex_name":cname,
                "model_num":model_num,
                "pdockq":if_stats['pdockq'],
                "ncontacts":if_stats['num_contacts'], 
                "plddt_min":round(if_stats['plddt'][0]),
                "plddt_avg":round(if_stats['plddt'][1]),
                "plddt_max":round(if_stats['plddt'][2]),
                "pae_min":round(if_stats['pae'][0]),
                "pae_avg":round(if_stats['pae'][1]),
                "pae_max":round(if_stats['pae'][2]), 
                "distance_avg":if_stats['distance_avg'],
            })


        stats = summarize_interface_statistics(interface_contacts)
        stats['best_model_num'] = best_interface_stats['model_num']
        stats['best_pdockq'] = best_interface_stats['pdockq']
        stats['best_plddt_avg'] = best_interface_stats['plddt'][1]
        stats['best_pae_avg'] = best_interface_stats['pae'][1]
        summary_stats[cname] = stats

        print("Finished analyzing " + cname)


    if len(summary_stats) < 1:
        print("Was not able to generate any summary statistics")
        return
    
    #output all the calculated values as CSV files into the specifed output folder (indexed by CPU to avoid different threads overwriting eachother)

    summary_df = pd.DataFrame.from_dict(summary_stats, 
                                        orient='index', 
                                        columns = ['avg_n_models',
                                                    'max_n_models', 
                                                    'num_contacts_with_max_n_models', 
                                                    'num_unique_contacts', 
                                                    'best_model_num', 
                                                    'best_pdockq',
                                                    'best_plddt_avg',
                                                    'best_pae_avg'])
    
    summary_df.index.name = 'complex_name'
    summary_df.to_csv(os.path.join(output_folder, f"summary_cpu{cpu_index}.csv"))

    interfaces_df = pd.DataFrame(all_interface_stats)
    interfaces_df.to_csv(os.path.join(output_folder, f"interfaces_cpu{cpu_index}.csv"), index=None)

    #there are cases when no contacts may be detected where we don't want to output anything
    if len(all_contacts) > 0:
        contacts_df = pd.DataFrame(all_contacts)
        contacts_df.to_csv(os.path.join(output_folder, f"contacts_cpu{cpu_index}.csv"), index=None)


def analysis_thread_did_finish(arg1):
   return


def analyze_folder(data_folder:str, max_distance:float, plddt_cutoff:float, pae_cutoff:float, pae_mode:str) -> str:
    """
        Analyze a folder containing protein structures in PDB format.

        :param data_folder (str): The path to the folder containing PDB and pAE JSON files to be analyzed.
        :param max_distance (float): The maximum distance between two atoms for them to be considered in contact.
        :param plddt_cutoff (float): The minimum pLDDT score required for a residue to be considered "well-modeled".
        :param pae_cutoff (float): The maximum predicted atomic error allowed for a residue to be considered "well-modeled".
        :param pae_mode (str): The method to use for calculating predicted atomic error (PAE). Possible values are "min" or "avg".
    """

    complex_names = get_finished_complexes(data_folder)
    if len(complex_names) < 1:
        print("ERROR: No complexes to analyze found. Please ensure all finished complexes/predictions you would like analyzed have a .done.txt file")
        return
    
    print(f"Found {len(complex_names)} complexes to analyze in folder: {data_folder}")

    output_folder = os.path.basename(data_folder) + "_analysis"
    index = 1
    while os.path.isdir(output_folder):
        #if we find existing folders with the output folder name we will iterate over index until we find an unused folder name
        output_folder = os.path.basename(data_folder) + "_analysis_" + str(index)
        index += 1

    #guaranteed to have a new unique output_folder name, lets make it     
    os.mkdir(output_folder)

    #find how many CPUs the system has and use as many as possible to speed up the analysis time
    num_cpus_to_use = mp.cpu_count()
    num_cpus_to_use = min(num_cpus_to_use,len(complex_names))
    print(f"Splitting analysis job across {num_cpus_to_use} different CPUs")
    pool = mp.Pool(num_cpus_to_use)

    #take the list of complexes and divide it across as many CPUs as we found
    complex_name_lists = distribute(complex_names, num_cpus_to_use)

    traces = []
    for cpu_index in range(0, num_cpus_to_use):
        #create a new thread to analyze the complexes (1 thread per CPU)
        traces.append(pool.apply_async(analyze_complexes, 
                                       args=(cpu_index, data_folder, output_folder, complex_name_lists[cpu_index], max_distance, plddt_cutoff, pae_cutoff, pae_mode), 
                                       callback=analysis_thread_did_finish))

    for t in traces:
        t.get()

    pool.close()
    pool.join()

    #merge all the seperate files produced by the independently running CPU threads
    for name in ['summary', 'interfaces', 'contacts']:

        files = glob.glob(os.path.join(output_folder, name + '_cpu*.csv'))
        if(len(files) < 1):
            # no files to join, skip
            continue

        sort_col = None
        if name == 'summary':
            #sort summary files by average num models descending as a default to bring top/most confident hits to the top
            sort_col = 'avg_n_models'
        join_csv_files(files, os.path.join(output_folder, name + '.csv'), sort_col=sort_col)

        #delete all the non-merged files produced by all the seperate CPU threads
        [os.remove(f) for f in files]

    return output_folder


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument(
        "input",
        default="",
        help="One or more folders with PDB files and pAE JSON files output by Colabfold. Note that '.done.txt' marker files produced by Colabfold are used to find the names of complexes to analyze.",
        nargs='*',)
    parser.add_argument(
        "--distance",
        default=8,
        help="Maximum distance in Angstroms that any two atoms in two residues in different chains can have for them be considered in contact for the analysis. Default is 8 Angstroms.",
        type=float,)
    parser.add_argument(
        "--pae",
        default=15,
        help="Maximum predicted Angstrom Error (pAE) value in Angstroms allowed for a contact(pair of residues) to be considered in the analysis. Valid values range from 0 (best) to 30 (worst). Default is 15.",
        type=float,)
    parser.add_argument(
        "--pae-mode",
        default='min',
        help=" How to combine the dual pAE values (x, y) and (y, x) outout for each residue pair in the pAE JSON files into a single pAE value for a residue pair (x, y).  Default is 'min'.",
        type=str, 
        choices=['min', 'avg'])
    parser.add_argument(
        "--plddt",
        default=50,
        help="Minimum pLDDT values required by both residues in a contact in order for that contact to be included in the analysis. Values range from 0 (worst) to 100 (best). Default is 50",
        type=float,)
    parser.add_argument(
        '--combine-all', 
        help="Combine the analysis from multiple folders specified by the input argument",
        action='store_true')
    args = parser.parse_args()

    if (args.distance < 1):
        sys.exit("The distance cutoff has been set too low. Please use a number greater than 1 Angstrom")

    if (args.pae < 1):
        sys.exit("The pAE cutoff has been set too low. Please use a number greater than 1 Angstrom")

    if (args.plddt < 1):
        sys.exit("The pLDDT cutoff has been set too low. Please use a number greater than 1")

    if (args.plddt > 99):
        sys.exit("The pLDDT cutoff has been set too high (pLDDT values range from 0 to 100). Please use a number less than 100 ")

    if len(args.input) < 1:
        sys.exit("No folders to analyze were provided");

    #loop through all the folders specified in the input
    output_folders = []
    for folder in args.input:
        
        if not os.path.isdir(folder):
            print(f"ERROR {folder} does not appear to be a non valid folder, skipping")
            continue

        print(f"Starting to analyze folder ({folder})")
        output_folders.append(analyze_folder(folder, args.distance, args.plddt, args.pae, args.pae_mode))
        print(f"Finished analyzing folder ({folder})")
        print(" "*80)
        print("*"*80)
        print("*"*80)
        print(" "*80)

    if len(args.input) > 1 and args.combine_all:

        combined_output_folder = 'combined_analysis_output'
        os.mkdir(combined_output_folder)

        for name in ['summary', 'interfaces', 'contacts']:
            
            csv_files = []
            for folder in output_folders:
                csv_files.append(os.path.join(folder, name + '.csv'))

            sort_col = None
            if name == 'summary':
                sort_col = 'avg_n_models'
            
            join_csv_files(csv_files, os.path.join(combined_output_folder, name + '.csv'), sort_col=sort_col)
        
        for folder in output_folders:
            shutil.rmtree(folder)

    print(f"Finished analyzing all specified folders")
