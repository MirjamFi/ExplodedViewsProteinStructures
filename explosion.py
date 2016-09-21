# run C:/Users/Figaschewski/Dropbox/Masterarbeit/Masterthesis/explosion.py

from pymol import cmd 
import get_colors
import math
import center_of_mass as cenma

cmd.delete('all')
cmd.set('ignore_case', 'off') ## case sensitive for chain ids

''' show axes '''
# import draw_axes as da
# da.axes()

''' load complex from PDB'''
pdbName = '3oaa'
pdbFile = pdbName + '.pdb1'
cmd.fetch(pdbName, type='pdb1')
cmd.set('all_states', 'on')
cmd.remove('solvent')
cmd.show('spheres', 'organic')

''' get chains of complex'''
chains = cmd.get_chains(pdbName)

#'''create object for each chain in selection '''
# cmd.split_chains(pdbName)

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

''' translate an object relative to complex of origin using center of mass'''
def translate_selection(original_xyz, trans_xyz, transname, factor):
	## vector between COMs to translate chain
	dist=math.sqrt((trans_xyz[0]-original_xyz[0])**2 +
					(trans_xyz[1]-original_xyz[1])**2 +
					(trans_xyz[2]-original_xyz[2])**2)

	vector=((trans_xyz[0]-original_xyz[0])/dist, 
			(trans_xyz[1]-original_xyz[1])/dist, 
			(trans_xyz[2]-original_xyz[2])/dist)
	trans_vec = [x * factor for x in vector] 

	## translate chain with vector
	cmd.translate(trans_vec,transname, camera=0)

''' Explosion'''
## create models for all states of complex
states = cmd.split_states(pdbName)

"""for each state create an object of every chain. In this state
 calulate COM for state and according chains and translate chain along 
 vector between COMs of state and chain"""
 
com_num = ""
for state in range(1,cmd.count_states()+1):
	## state = 0 would calculate the COM for all states at once
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
		
		## only translate chain if there are multiple chains
		if len(chains) > 1:	
			## translate chain
			translate_selection(complex_xyz, chain_xyz, chainname, 50)
			cmd.delete(chainname + '_COM')
		
		## translate ligands on current chain, if ligands are present
		if stored.ligands:
		
			## if chain has been translated, calculate COM of translated chain
			if len(chains) > 1:
				cenma.com(chainname, state = 1)
				cmd.zoom(chainname + '_COM')
				chain_xyz = cmd.get_position(chainname + '_COM')
					
			## delete chain's COM
			cmd.delete(chainname + '_COM')
			
			## select ligand on chain
			selection = chainname + " & organic"
			ligandname = chainname + "_org"
			
			## only translate if chain really has a ligand
			if cmd.count_atoms(selection) > 0:
				cmd.create(ligandname, selection)
				
				## calculate COM of ligand
				cenma.com(ligandname, state = 1)
				cmd.zoom(ligandname + '_COM')
				ligand_xyz = cmd.get_position(ligandname + '_COM')
				
				## translate ligand
				translate_selection(chain_xyz, ligand_xyz, ligandname, 30)

				## delete ligand's COM
				cmd.delete(ligandname + '_COM')
			
				## color binding site
				binding = 'byres (' + chainname + ' nto. 3 of organic)'
				cmd.select('inter', binding)
				cmd.color('red', 'inter')
				cmd.delete('inter')
				
				## remove ligand from translated chain
				cmd.remove(selection)	
				
		## delete chain's COM
		cmd.delete(chainname + '_COM')
	
	## delete state object
	cmd.delete(pdbName + "_"+ s)
	
	## delete state's COM
	cmd.delete(pdbName + "_" + s + '_COM')
		
cmd.zoom('all')
			