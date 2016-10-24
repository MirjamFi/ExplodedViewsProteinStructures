# Exploded Views for Molecular Structures

## Abstract

3D protein structure models such as those deposited in the protein databank (PDB) often consist of multiple parts, e.g. multiple proteins, chains, binding partners (DNA, RNA), and ligands. Sometimes, important parts are occluded from the default viewpoint or even from *any* camera position. In order get an overview of the composition of such complex arrangements, *exploded views* as known from engineering will be applied to molecular structures in this work. An illustrative example of such an exploded view for the 70S Ribosome that was created manually is shown here on the right:
<img src="https://bmcresearch.utm.utoronto.ca/sciencevislab/wp-content/media/2013/12/70S_Ribosome_1920x1080.jpg" style="width: 640px;"></img>

We anticipate that interactive exploded views will greatly help users (in particular students) to (1) get an overview of composition of complex molecular structures and (2) better understand how the individual parts of a complex work. An additional benefit of this work is to quickly to obtain overview figures that can be used for presentation.

## Related Work
[Bruckner and Groeller](http://dl.acm.org/citation.cfm?id=1187828) first describe an algorithm for the computation of exploded views for volumetric data in a medical context. [Li et al.](http://dl.acm.org/citation.cfm?id=1360700) present methods for the automatic generation of exploded view diagrams for geometric objects in an engineering context. 

## Research Questions / Contributions
* Which are the *canonical explosion directions* for molecule complexes?
* How to automatically label all parts without occlusion?
* How to determine *parts* and their hierarchy of a molecular complex? 
* How to determine the *order* by which parts are exploded?
* How to interactively steer this process?
* How to handle dynamic parts?

## Desired Features
* User defines a rendering mode, e.g. cartoon, surface, or ball and stick
* User selects a 'part' to expose (e.g. ligand or protein in complex, center of mass)
* explode interactively, i.e. let user control speed
* Label parts if not exploding
* Select good viewpoint(s) for each frame
* Export to movie

## Project Outline

1. Find and read related work
2. Become familiar with test structures, PDB files, and PyMOL
3. Sketch exploding algorithm for simple test case, e.g. [Chaperonin](http://www.rcsb.org/pdb/explore.do?structureId=1AON)
4. Implement prototype in Python/PyMOL
5. Convey user study (e.g. with Biochemistry students) to see if exploded views support understanding of complex structures
6. Write thesis

## Implementation Steps
1. Define part hierarchy and how to obtain it from PDB files
	* Complex, Ligands, Biological Unit, Asymetric Unit, Chains, SS, Atoms
2. For every part: Determine explosion direction
3. Compute explosion graph to determine order of parts
4. Resolve interlocks, if any
5. Explode (and compute best view)
6. Add labels

## Notes

### Li2008: 
* cutaways don't apply to Molecules
* molecules have natural *part hierarchies*
* [Illustrative Example](https://bmcresearch.utm.utoronto.ca/sciencevislab/wp-content/media/2013/12/70S_Ribosome_1920x1080.jpg)

## Example Structures that could be used for testing
* [Chaperonin](http://www.rcsb.org/pdb/explore.do?structureId=1AON)
* [2XNV](http://www.rcsb.org/pdb/explore/explore.do?structureId=2xnv)
* [3OAA](http://www.rcsb.org/pdb/explore/explore.do?structureId=3OAA),
* [Dengue Virus](http://www.rcsb.org/pdb/explore/explore.do?structureId=1k4r) (use fetch 1k4r, type=pdb1 to get biological unit and set all_states, on)
* PMP structures, molecular injection complex (e.g. http://www.ebi.ac.uk/pdbe/entry/emdb/EMD-6330)
* Splicosome, Proteasome, Polyokapsid, http://pdb101.rcsb.org/motm

Charlotta:
3oaa: E.coli ATP synthase (8-mer)
2IH3: K+ channel (12-mer)
5K7L: ein anderer K+ channel (8-mer)
1VQP: Transition state analogue "RAP" bound to the large ribosomal subunit (27-mer)
5EN5: bacterial efflux pump (6-mer)
5J4B: bacterial ribosome bound to cisplatin (27-mer)
5IT8: bacterial ribosome (55-mer)
1ERI: DNA + endonuclease (2-mer + DNA)
4ESV: DNA + helicase (6-mer)
3OM3: Tobacco Mosaic Virus (49-mer)

Für die ganz Ambitionierten:
2BUK: Satellite Tobacco Necrosis Virus (60-mer)
3IYN: cryo-EM structure of human adenovirus (1200-mer)
1OHG: capsid (420-mer)

4V59: fatty acid synthase complex

Jens:
Pore in Membran

## Applications
### Labeling
This is from a [paper](http://www.biodiscoveryjournal.co.uk/Archive/A27.htm):
<img src='http://www.biodiscoveryjournal.co.uk/Archive/Media/A27Figure-5.jpg'></img>

