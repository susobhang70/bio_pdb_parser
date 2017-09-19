# PDB Parser

Script that takes pdbfile name as input to retrieve the following information and write to 2wsc_output.txt.

- Title/Name  of the protein.
- Total length of the protein/ number of residues in the given pdb 
- No of chains present in protein and their names(Ascending order if chains are named numerically followed by alphabetical order)
- All aminoacid ratios present in the protein in alphabetical order. Ex: Ratio(Leu) = no of leucine present/total length of protein.
- If there any unknown aminoacids present and if so, mentioned total count.
- If there any ligand molecules other than water and if so, mentioned their names.
- Calculated all possible phi, psi, omega angles for the given pdb.
