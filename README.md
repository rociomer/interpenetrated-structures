# interpenetrated-structures

## Description
This program was written so as to facilitate the generation 
of interpenetrated structures (2-fold and greater) of covalent 
organic frameworks (COFs) and similar nanoporous frameworks.

## Instructions
(Note: will need Zeo++ to convert files from CSSR to CIF format)

Copy *file_conversion.py* and *find_interp_structures.py* to
directory with CSSR files for structures to interpenetrate. If necessary,
modify the path to Zeo++ at the top of *file_conversion.py*, and modify any 
parameters at the top of *find_interp_structures.py* (such as 
the displacement length and overlap distance between the interpenetrated
structures to be generated). Then simply run *find_interp_structures.py* 
(uses Python2.7). Output files will be written to the working directory.

For sample cluster submission scripts, see *cluster-submission-scripts/*.

## Authors
Rocío Mercado

Ray Fu

## Link 
https://github.com/rociomer/interpenetrated-structures

## Licence
None

