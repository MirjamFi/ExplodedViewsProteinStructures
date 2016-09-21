# run C:/Users/Figaschewski/Dropbox/Masterarbeit/test.py

from pymol import cmd 
#from Bio.PDB import * # parse PDB file
#import elements_from_PDB as efPDB 
import get_colors
import math
import center_of_mass as cenma

cmd.delete('all')
cmd.set('ignore_case', 'off') ## case sensitive for chain ids

''' show axes '''
# import draw_axes as da
# da.axes()

''' load complex from PDB'''
pdbName = '1k4r'
pdbFile = pdbName + '.pdb1'
cmd.fetch(pdbName, type='pdb1')
cmd.set('all_states', 'on')
cmd.remove('solvent')

''' get chains of complex'''
chains = cmd.get_chains(pdbName)

'''create object for each chain in selection '''
# cmd.split_chains(pdbName)

''' parse pdb '''
# parser = PDBParser()
# structure = parser.get_structure(pdbName, pdbFile)
# chains = []

''' get all ligands (organic) from complex'''
"""import stored for passing data back and forth"""
from pymol import stored
stored.ligands = []

## create the original selection of all organic atoms 
cmd.select("allOrg", "org")

##" temporary selection
cmd.select("temp", "none")

""" while initial selection is not empty, "pop" from it and query the atom's 
 molecule's name """
while cmd.count_atoms("allOrg") != 0:

	## pop--peek, rather
	cmd.select("temp", "first allOrg")
	
	## store the residue name
	cmd.iterate("temp", "stored.ligands.append(resn)")
	
	## remove by molecule
	cmd.select("allOrg", "allOrg and not (bm. temp)")
	
cmd.delete('allOrg')
cmd.delete('temp')


''' Explosion'''
## create models for all states of complex
states = cmd.split_states(pdbName)

"""for each state create an object of every chain. In this state
 calulate COM for state and according chains and translate chain along 
 vector between COMs of state and chain"""
com_num = ""
for state in range(1,cmd.count_states()+1):
	if state == 0:
		continue

	## name chain object in this state
	s = "00" + str(state)
	if state < 10:
		s = "0" + s
		
	## COM of state (all chains)
	cenma.com(pdbName + "_" + s, state = 1)
	
	## get coordinates of COM
	cmd.zoom(pdbName + "_" + s + '_COM')
	complex_xyz = cmd.get_position(pdbName + '_COM' + com_num)
	
	## calculate COM for single chains in state
	for c in chains:
	
		## color each chain individually
		cmd.color(get_colors.get_random_color(), pdbName + "_" + s + "& chain " + c)

		## create object for chain with name
		chainname = s + '_'+c

		origin = pdbName+'_'+s+ ' & chain '+c
		cmd.create(chainname, origin)

		## COM of chain
		cenma.com(chainname, state = 1)
		cmd.zoom(chainname + '_COM')
		chain_xyz = cmd.get_position(chainname + '_COM')
		
		## vector between COMs to translate chain
		dist=math.sqrt((chain_xyz[0]-complex_xyz[0])**2 +
						(chain_xyz[1]-complex_xyz[1])**2 +
						(chain_xyz[2]-complex_xyz[2])**2)

		vector=((chain_xyz[0]-complex_xyz[0])/dist, 
				(chain_xyz[1]-complex_xyz[1])/dist, 
				(chain_xyz[2]-complex_xyz[2])/dist)
		trans_vec = [x * 50 for x in vector] 

		## translate chain with vector
		cmd.translate(trans_vec,chainname, camera=0)
		
		## delete chain's COM
		cmd.delete(chainname + '_COM')
		
	## delete state object
	cmd.delete(pdbName + "_"+ s)
	
	## delete state's COM
	cmd.delete(pdbName + "_" + s + '_COM')
		
cmd.zoom('all')



''' get all elements in pdb (molecule names with chains, chains with ligands) '''
# chains_per_molecule, ligands_per_chain = efPDB.getElementsFromPDB(structure,
										# pdbFile)	
# print chains_per_molecule, ligands_per_chain

# ''' show molecules in different colors and highlight ligands as spheres'''
# for key in chains_per_molecule:
	# if chains_per_molecule[key]:
		# chains = "+".join(chains_per_molecule[key])
		# chains = "chain " + chains
		# cmd.create(key, chains)
		# cmd.color(get_colors.get_random_color(), chains)
		
# cmd.show('spheres', 'organic')
# #cmd.delete(pdbName)	# delete original pdb, keep separated molecules

# ''' color binding site of ligands and coresponding chain'''
# for key in ligands_per_chain:
	# if ligands_per_chain[key]:
		# for value in ligands_per_chain[key]:
			# selection = 'byres (chain ' + key + ' nto. 3 of resn ' + value + ')'
			# cmd.select('inter', selection)
			# cmd.color('red', 'inter')
			# cmd.delete('inter')
			