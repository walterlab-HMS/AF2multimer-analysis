# AF2multimer-analysis


This python script allows one to find contacts between residues in multimeric structure files produced as output from Alphafold2 via the Colabfold pipeline <https://github.com/sokrypton/ColabFold/tree/main/colabfold>. It integrates both physical proximity and Alphafold conficence metrics such as the predicted Alignment Error(pAE) and the predicted Local Distance Difference Test (pLDDT) to determine whether a pair of residues is a valid contact.

## Running

```
usage: colabfold_analysis.py [-h] [--distance DISTANCE] [--pae PAE] [--pae-mode {min,avg}]
                             [--plddt PLDDT] [--combine-all]
                             [input [input ...]]

positional arguments:
  input                 One or more folders with PDB files and pAE JSON files output by Colabfold.
                        Note that '.done.txt' marker files produced by Colabfold are used to find
                        the names of complexes to analyze.

optional arguments:
  -h, --help            show this help message and exit
  
  --distance DISTANCE   Maximum distance in Angstroms that any two atoms in two residues in
                        different chains can have for them be considered in contact for the
                        analysis. Default is 8 Angstroms.
                        
  --pae PAE             Maximum predicted Angstrom Error (pAE) value in Angstroms allowed for a
                        contact(pair of residues) to be considered in the analysis. Valid values
                        range from 0 (best) to 30 (worst). Default is 15.
                        
  --pae-mode {min,avg}  How to combine the dual PAE values (x, y) and (y, x) into a single PAE
                        value for a residue pair (x, y). Default is 'min'.
                        
  --plddt PLDDT         Minimum pLDDT values required by both residues in a contact in order for
                        that contact to be included in the analysis. Values range from 0 (worst)
                        to 100 (best). Default is 50.
                        
  --aas AAS             A string representing what amino acids contacts to look/filter for. Allows you
                        to limit what contacts to include in the analysis. By default is blank meaning
                        all amino acids. A value of K would be for any lysine lysine pairs. KR would be
                        RR, KR, RK, or RR pairs, etc
                        
  --name-filter NAME_FILTER
                        An optional string that allows one to only analyze complexes that contain
                        that string in their name
                        
  --combine-all         Combine the analysis from multiple folders specified by the input argument
  
  --ignore-pae          Ignore PAE values and just analyze the PDB files. Overides any other PAE
                        settings.
```


Examples:

```
python3 colabfold_analysis.py my_exciting_colabfold_output_folder

python3 colabfold_analysis.py my_exciting_colabfold_output_folder --pae 12 --plddt 50 --pae-mode avg

python3 colabfold_analysis.py folder1 folder2 folder3 --pae 12 --plddt 50 --pae-mode avg --combine-all

python3 colabfold_analysis.py folder1 --aas DEHKR

python3 colabfold_analysis.py folder1 --ignore-pae --name-filter MCM

python3 colabfold_analysis.py folder_? --distance 10 --plddt 60 --pae-mode min --combine-all

```

## Example folder to be analyzed
Note that presence of a .done.txt marker file produced by Colabfold. Only complexes with a marker file will be analyzed.
```
colabfold_output/
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN.done.txt
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa.a3m
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa.done.txt
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_coverage.png
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_scores_rank_001_alphafold2_multimer_v3_model_5_seed_000.json
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_scores_rank_002_alphafold2_multimer_v3_model_4_seed_000.json
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_scores_rank_003_alphafold2_multimer_v3_model_2_seed_000.json
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_scores_rank_004_alphafold2_multimer_v3_model_1_seed_000.json
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_scores_rank_005_alphafold2_multimer_v3_model_3_seed_000.json
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_template_domain_names.json
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000.pdb
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_unrelaxed_rank_002_alphafold2_multimer_v3_model_4_seed_000.pdb
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_unrelaxed_rank_003_alphafold2_multimer_v3_model_2_seed_000.pdb
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000.pdb
├── MCM2_HUMAN__MCM5_HUMAN-CDC45_HUMAN__2204aa_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000.pdb
├── cite.bibtex
├── config.json
└── log.txt
```

## Outputs

Running this script will produce one or more folders each containing 3 comma seperated value (CSV) files that you can then open with a standard text editor or any spreadhseet program. 

The 3 files are: summary.csv, interfaces.csv, and contacts.csv. 

### summary.csv

Summarizes all the findings per complex across all models that were run for it. 
Each row is a summary for one complex.

| complex_name         | avg_n_models                         | max_n_models                                     | num_contacts_with_max_n_models                                         | num_unique_contacts                                      | best_model_num                                                                   | best_pdockq                                                                   | best_plddt_avg                                                                             | best_pae_avg                                                                             |
|----------------------|--------------------------------------|--------------------------------------------------|------------------------------------------------------------------------|----------------------------------------------------------|----------------------------------------------------------------------------------|-------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| name of the complex  | avg number of <br>models per contact | max number of models any <br>contact was seen in | number of unique contacts that<br> were seen max model number of times | number of unique contacts <br>across all models anlayzed | model number of prediction<br>producing strongest interaction <br>score (pdockq) | highest pdockq score recorded <br>across all predictions for this <br>complex | the average pLDDT values <br>across the interface for the<br>model with the highest pDOCKQ | the average pAE values <br>across the interface for the<br>model with the highest pDOCKQ |

### interfaces.csv
Shows the statistics for each prediction made for each complex. 
Each row is 1 prediction (structure/JSON score file)

| complex_name         | model_num       | pdockq                                                                    | ncontacts                             | plddt_min                                        | plddt_avg                                    | plddt_max                                        | pae_min                                        | pae_avg                                            | pae_max                                        | distance_avg                                                             |
|----------------------|-----------------|---------------------------------------------------------------------------|---------------------------------------|--------------------------------------------------|----------------------------------------------|--------------------------------------------------|------------------------------------------------|----------------------------------------------------|------------------------------------------------|--------------------------------------------------------------------------|
| name of the complex  | AF model number | predicted DOCKQ interface accuracy score<br>ranges from 0 worst to best 1 | number of contacts seen in prediction | Min residue pair pLDDT observed in the interface | Average pair pLDDT observed in the interface | Max residue pair pLDDT observed in the interface | Min residue pair PAE observed in the interface | Average residue pair PAE observed in the interface | Max residue pair PAE observed in the interface | Average distance between closest atoms in residue pairs in the interface |


### contacts.csv
A comprehensive table of all residue contact pairs between all chains that met the contact criteria specified during the run. 
Each row is 1 pair of interacting residues in different chains.

| complex_name        | model_num              | aa1_chain             | aa1_index                           | aa2_chain             | aa1_plddt     | aa2_index                           | aa2_type                    | aa2_plddt     | aa1_type                    | pae                                                                       | min_distance                                          |
|---------------------|------------------------|-----------------------|-------------------------------------|-----------------------|---------------|-------------------------------------|-----------------------------|---------------|-----------------------------|---------------------------------------------------------------------------|-------------------------------------------------------|
| Name of the complex | AlphaFold model number | chain residue 1 is in | Index of residue 1 within its chain | chain residue 2 is in | pLDDT for aa1 | Index of residue 2 within its chain | 1 letter code for residue 2 | pLDDT for aa2 | 1 letter code for residue 1 | Combined pAE value for residue pair calculated using specified "pae_mode" | Minimum distance in angstroms between the 2 residues. |


