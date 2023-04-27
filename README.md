# AF2multimer-analysis
multimer-analysis


This python script allows one to find contacts between residues in multimeric structure files produced as output from Alphafold2 via the Colabfold pipeline. 

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
  --distance DISTANCE   Maximum distance in Angstroms for a contact to be included in the final
                        output. Default is 8. 
  --pae PAE             Maximum PAE value allowed for a contact to be included in the final output. Valid values from 0 to 30. Default is 15.
  --pae-mode {min,avg}  How to combine the dual PAE values (x, y) and (y, x) single PAE value for
                        a residue pair (x, y) Valid values from 0 to 30.  Default is "min"
  --plddt PLDDT         Minimum pLDDT values required by both residues in a contact in order for
                        that contact to be included in the final output. Valid values from 0 to 100. Default is 70.
  --combine-all         Combine the analysis from multiple folders specified by the input argument
```
