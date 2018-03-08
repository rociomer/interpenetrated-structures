# interpenetrated-structures

## Authors
Rocío Mercado

Ray Fu

## Link 
https://github.com/rociomer/interpenetrated-structures

## Licence
None

## Instructions
(Note: will need Zeo++ to convert files from CSSR to CIF format)

Copy fileConversion.py and findInterpenetratedStructures.py to
directory with CSSR files for structures to interpenetrate. If necessary,
modify the path to Zeo++ at the top of fileConversion.py, and modify any 
parameters at the top of findInterpenetratedStructures.py (such as 
the displacement length and overlap distance between the interpenetrated
structures to be generated). Then simply run findInterpenetratedStructures.py 
(uses Python2.7). Output files will be written to the working directory.

For sample cluster submission scripts, see *cluster-submission-scripts/*.

## Description
This program was written so as to facilitate the generation 
of interpenetrated structures (2-fold and greater) of covalent 
organic frameworks (COFs) and similar nanoporous frameworks.
