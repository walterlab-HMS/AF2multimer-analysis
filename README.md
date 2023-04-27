# AF2multimer-analysis


This python script allows one to find contacts between residues in multimeric structure files produced as output from Alphafold2 via the Colabfold pipeline <https://github.com/sokrypton/ColabFold/tree/main/colabfold>. It integrates both physical proximity and Alphafold conficence metrics such as the predicted Alignment Error(pAE) and the predicted Local Distance Difference Test (pLDDT) to determine whether residues are for a valid contact.


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
                        
  --combine-all         Combine the analysis from multiple folders specified by the input argument
```


Examples to run:

```
python3 colabfold_analysis.py my_exciting_colabfold_output_folder

python3 colabfold_analysis.py my_exciting_colabfold_output_folder --pae 12 --plddt 50 --pae-mode avg

python3 colabfold_analysis.py folder1 folder2 folder3 --pae 12 --plddt 50 --pae-mode avg --combine-all

python3 colabfold_analysis.py folder_? --distance 10 --plddt 60 --pae-mode min --combine-all

```
