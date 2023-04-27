# AF2multimer-analysis
multimer-analysis


This python script allows one to find contacts between residues in multimeric structure files produced as output from Alphafold2 via the Colabfold pipeline. <https://github.com/sokrypton/ColabFold/tree/main/colabfold>

```
usage: colabfold_analysis.py [-h] [--distance DISTANCE] [--pae PAE] [--pae-mode {min,avg}]
                             [--plddt PLDDT] [--combine-all]
                             [input [input ...]]

positional arguments:
  input                 One or more folders with PDB files and pAE json files output by Colabfold.
                        Note the '.done.txt' files are used to find the names of complexes to
                        analyze.

optional arguments:
  -h, --help            show this help message and exit
  
  --distance DISTANCE   Maximum distance in Angstroms that any two atoms in two residues in
                        different chains can have for them be considered in contact for the
                        analyis. Default is 8 Angstroms.
                        
  --pae PAE             Maximum predicted Angstrom Error (pAE) value in Angstroms allowed for a
                        contact(pair of residues) to be considered in the analysis. Valid values
                        range from 0 (best) to 30 (worst). Default is 15.
                        
  --pae-mode {min,avg}  How to combine the dual PAE values (x, y) and (y, x) into a single PAE
                        value for a residue pair (x, y). Default is 'min'.
                        
  --plddt PLDDT         Minimum pLDDT values required by both residues in a contact in order for
                        that contact to be included in the analysis.
                        
  --combine-all         Combine the analysis from multiple folders specified by the input argument
```
